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


precice_root = os.getenv('PRECICE_ROOT')
precice_python_adapter_root = precice_root+"/src/precice/adapters/python"
sys.path.insert(0, precice_python_adapter_root)

import PySolverInterface
from PySolverInterface import *

print "Starting Left Solver..."

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)
args = parser.parse_args()

style = config.style
configFileName = args.configurationFileName
N = config.n_elem
init_temp = config.init_temp
bound_temp = config.bound_temp
jump_temp = config.jump_temp
dt = config.dt
maxT = max(init_temp, bound_temp, jump_temp)

solverName = "TEMPLEFT"

print "Configure preCICE..."
interface = PySolverInterface(solverName, 0, 1)
interface.configure(configFileName)

dimensions = interface.getDimensions()

leftNodes = (init_temp * np.ones(N*N)).reshape(N,N)
leftNodes = np.pad(leftNodes, (1,1), 'constant', constant_values=(bound_temp,bound_temp))[0:N+1][:]
leftNodes_n = np.copy(leftNodes)

midNodes = np.copy(leftNodes[N][:])
rightBound = np.copy(midNodes)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

meshID = interface.getMeshID("leftNodes")
leftDataID = interface.getDataID("tempMid", meshID)
rightDataID = interface.getDataID("rightBound", meshID)

vertexIDs = np.zeros(N+2)
grid = np.zeros([dimensions, N+2])

grid[0,:] = np.linspace(0, 1, N+2)

interface.setMeshVertices(meshID, N+2, grid.flatten('F'), vertexIDs)

t = 0

precice_tau = interface.initialize()

if interface.isActionRequired(PyActionWriteInitialData()):
    interface.writeBlockScalarData(leftDataID, N+2, vertexIDs, tempMid)
    interface.fulfilledAction(PyActionWriteInitialData())

interface.initializeData()

if interface.isReadDataAvailable():
    interface.readBlockScalarData(rightDataID, N+2, vertexIDs, rightBound)

while interface.isCouplingOngoing():
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())
    
    if style == "explicit":
    	leftNodes_n, midNodes = temp_step_ex(leftNodes, midNodes, rightBound, N, dt, 1)
    else:
    	leftNodes_n, midNodes = temp_step(leftNodes, midNodes, rightBound, N, dt, 1)

    interface.writeBlockScalarData(leftDataID, N+2, vertexIDs, midNodes)
    interface.advance(precice_tau)
    interface.readBlockScalarData(rightDataID, N+2, vertexIDs, rightBound)

    if interface.isActionRequired(PyActionReadIterationCheckpoint()): # i.e. not yet converged
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else:
        t += dt
        tubePlotting.doPlottingLeft(ax, midNodes, leftNodes, t, maxT)
        ax.cla()
        leftNodes = np.copy(leftNodes_n)

print "Exiting Left Solver"

interface.finalize()

