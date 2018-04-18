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

print "Starting Right Solver..."

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)
args = parser.parse_args()

style = config.style
configFileName = args.configurationFileName
N = config.n_elem
init_temp = config.init_temp
bound_temp = config.bound_temp
jump_temp =config.jump_temp
dt = config.dt
maxT = max(init_temp, bound_temp, jump_temp)

solverName = "TEMPRIGHT"

print "Configure preCICE..."
interface = PySolverInterface(solverName, 0, 1)
interface.configure(configFileName)

dimensions = interface.getDimensions()

rightNodes = (init_temp * np.ones((N-1)*N)).reshape((N-1),N)
rightNodes = np.pad(rightNodes, (1,1), 'constant', constant_values=(bound_temp,bound_temp))[0:N][:]
rightNodes_n = np.copy(rightNodes)

midNodes = np.copy(rightNodes[N-1][:])
rightBound = np.copy(midNodes)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

meshID = interface.getMeshID("rightNodes")
leftDataID = interface.getDataID("tempMid", meshID)
rightDataID = interface.getDataID("rightBound", meshID)

vertexIDs = np.zeros(N+2)
grid = np.zeros([dimensions, N+2])

grid[0,:] = np.linspace(0, 1, N+2)

interface.setMeshVertices(meshID, N+2, grid.flatten('F'), vertexIDs)

t = 0

precice_tau = interface.initialize()

if (interface.isActionRequired(PyActionWriteInitialData())):
   interface.writeBlockScalarData(rightDataID, N+2, vertexIDs, rightBound)
   interface.fulfilledAction(PyActionWriteInitialData())

interface.initializeData()

if (interface.isReadDataAvailable()):
   interface.readBlockScalarData(leftDataID, N+2, vertexIDs, midNodes)

while interface.isCouplingOngoing():
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())

    if style == "explicit":
      rightNodes_n, rightBound = temp_step_ex(rightNodes, rightBound, midNodes, N, dt, 0)
    else:
      rightNodes_n, rightBound = temp_step(rightNodes, rightBound, midNodes, N, dt, 0)

    interface.writeBlockScalarData(rightDataID, N+2, vertexIDs, rightBound)
    interface.advance(precice_tau)
    interface.readBlockScalarData(leftDataID, N+2, vertexIDs, midNodes)

    if interface.isActionRequired(PyActionReadIterationCheckpoint()):
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else:
        t += dt
        tubePlotting.doPlottingRight(ax, rightBound, rightNodes, t, maxT)
        ax.cla()
        rightNodes = np.copy(rightNodes_n)

print "Exiting Right Solver"

interface.finalize()

