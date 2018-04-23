from __future__ import division
from thetaScheme import temp_step
from thetaScheme import temp_step_ex
import os
import sys
import argparse
import tubePlotting
from mpi4py import MPI
import numpy as np
import configuration_file as config
import matplotlib.pyplot as plt

# setup of precice environment. locates the latest precice build and python adapter
precice_root = os.getenv('PRECICE_ROOT')
precice_python_adapter_root = precice_root+"/src/precice/adapters/python"
sys.path.insert(0, precice_python_adapter_root)

import PySolverInterface
from PySolverInterface import *

print "Starting Left Solver..."

# reads in precice configuration file
parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)
args = parser.parse_args()

# note all configuration variables are taken from configuration_file.py
style = config.style
N = config.n_elem
init_temp = config.init_temp
bound_temp = config.bound_temp
dt = config.dt
maxT = max(init_temp, bound_temp, jump_temp)

solverName = "TEMPLEFT"

print "Configure preCICE..."
# initialize the PySovlverInterface class. Only one instance is required as this is sequential
interface = PySolverInterface(solverName, 0, 1)
interface.configure(configFileName)

# reads dimensions from precice configuration file
dimensions = interface.getDimensions()

# initialization of local data
leftNodes = (init_temp * np.ones(N*N)).reshape(N,N)
leftNodes = np.pad(leftNodes, (1,1), 'constant', constant_values=(bound_temp,bound_temp))[0:N+1][:]
leftNodes_n = np.copy(leftNodes)

midNodes = np.copy(leftNodes[N][:])
rightBound = np.copy(midNodes)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# pull the mesh and data IDs used to read/write to/from precice
meshID = interface.getMeshID("leftNodes")
leftDataID = interface.getDataID("tempMid", meshID)
rightDataID = interface.getDataID("rightBound", meshID)

vertexIDs = np.zeros(N+2)
grid = np.zeros([dimensions, N+2])

grid[0,:] = np.linspace(0, 1, N+2)

# mesh geometry defined
interface.setMeshVertices(meshID, N+2, grid.flatten('F'), vertexIDs)

t = 0

# dt is taken from precice, as defined in precice configuration file
precice_tau = interface.initialize()

# data is initialized before start of simulation
if interface.isActionRequired(PyActionWriteInitialData()):
    interface.writeBlockScalarData(leftDataID, N+2, vertexIDs, tempMid)
    interface.fulfilledAction(PyActionWriteInitialData())

interface.initializeData()

if interface.isReadDataAvailable():
    interface.readBlockScalarData(rightDataID, N+2, vertexIDs, rightBound)
# while loop is used for simulation process. Auxilliary method isCouplingOngoing is used as a break condition
while interface.isCouplingOngoing():
    # basic structure here to allow easy implementation of checkpointing here. Checkpointing is inherent in the update
    # step here and does not use the preCICE API
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())
    
    # will perform explicit time-stepping (temp_step_ex) if style is set to explicit in configuration file
    # otherwise will perform implicit time-stepping. Remove or comment out if-statement to get significant performance
    # improvement if time-stepping method required is known.
    if style == "explicit":
    	leftNodes_n, midNodes = temp_step_ex(leftNodes, midNodes, rightBound, N, dt, 1)
    else:
    	leftNodes_n, midNodes = temp_step(leftNodes, midNodes, rightBound, N, dt, 1)

    # Updated data is written to preCICE. Unique VertexIDs are referenced to write data to correct location.
    interface.writeBlockScalarData(leftDataID, N+2, vertexIDs, midNodes)
    # Time step is advanced.
    interface.advance(precice_tau)
    interface.readBlockScalarData(rightDataID, N+2, vertexIDs, rightBound)

    # Basic structure here to allow easy implementation of checkpointing here. Checkpointing is inherent in the update
    # step here and does not use the preCICE API
    if interface.isActionRequired(PyActionReadIterationCheckpoint()): # i.e. not yet converged
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else:
        t += dt
        # plotting of data for each time step using matplotlib. comment out this line for significant performance increase
        # during run-time
        tubePlotting.doPlottingLeft(ax, midNodes, leftNodes, t, maxT)
        ax.cla()
        # changes to local data are preserved upon satisfying convergence requirements.
        leftNodes = np.copy(leftNodes_n)

print "Exiting Left Solver"


# finalize() called to perform cleanup and close communication channels
interface.finalize()

