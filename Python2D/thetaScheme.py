from __future__ import division
import numpy as np

import configuration_file as config
import math

# implicit time-stepping via Euler's
def temp_step(temp, nearBound, farBound, N, dt, solverNum):
    scaling = config.scalingFactor
    temp_n = np.copy(temp)

    nearBound_n = np.copy(nearBound)
    offDiag = -dt*(N*N)*scaling
    diag = 1 + 4*dt*(N*N)*scaling
    precision = 1e-2
    error = 10

    while(error > precision):
        last_temp = np.copy(temp_n)
        last_bound = np.copy(nearBound_n)
        for col in range(1, N+1):
            temp_n[-1][col] = (temp[-1][col] - offDiag*nearBound_n[col] - offDiag*temp_n[-2][col] 
                                - offDiag*temp_n[-1][col+1] - offDiag*temp_n[-1][col-1])/diag
            nearBound_n[col] = (nearBound[col] - offDiag*temp_n[-1][col] - offDiag*farBound[col] 
                                - offDiag*nearBound_n[col+1] - offDiag*nearBound_n[col-1])/diag

        for row in range(1, N+solverNum-1):
            for col in range(1, N+1):
                temp_n[row][col] = (temp[row][col] - offDiag*temp_n[row+1][col] - offDiag*temp_n[row-1][col] 
                                - offDiag*temp_n[row][col+1] - offDiag*temp_n[row][col-1])/diag
        
        error = math.sqrt(np.sum(np.power(temp_n - last_temp, 2)) + np.sum(np.power(nearBound_n - last_bound, 2)))
        print(error)

    return temp_n, nearBound_n

# explicit time-stepping via Euler's
def temp_step_ex(temp, nearBound, farBound, N, dt, solverNum):
    scaling = config.scalingFactor
    grad = np.zeros((temp.shape))
    grad_bound = np.zeros(nearBound.shape)
    diag = -4*N*N*dt*scaling
    offDiag = 1*N*N*dt*scaling

    for col in range(1, N+1):
        grad[-1][col] = diag*temp[-1][col] + offDiag*nearBound[col] + offDiag*temp[-2][col] + offDiag*temp[-1][col+1] + offDiag*temp[-1][col-1]
        grad_bound[col] = diag*nearBound[col] + offDiag*temp[-1][col] + offDiag*farBound[col] + offDiag*nearBound[col+1] + offDiag*nearBound[col-1]

    for row in range(1, N+solverNum-1):
        for col in range(1, N+1):
            grad[row][col] = diag*temp[row][col] + offDiag*temp[row+1][col] + offDiag*temp[row-1][col] + offDiag*temp[row][col+1] + offDiag*temp[row][col-1]
    temp = temp + grad
    nearBound = nearBound + grad_bound
    
    return temp, nearBound

def temp_step_whole(temp, N, dt):
    scaling = config.scalingFactor
    temp_n = np.copy(temp)

    offDiag = -dt*(N*N)*scaling
    diag = 1 + 4*dt*(N*N)*scaling
    precision = 1e-2
    error = 10

    while(error > precision):

        last_temp = np.copy(temp_n)
        for row in range(1, 2*N+2):
            for col in range(1, N+1):
                temp_n[row][col] = (temp[row][col] - offDiag*temp_n[row+1][col] - offDiag*temp_n[row-1][col] 
                                - offDiag*temp_n[row][col+1] - offDiag*temp_n[row][col-1])/diag
        error = math.sqrt(np.sum(np.power(temp_n - last_temp, 2)))

    return temp_n

def temp_step_whole_ex(temp, N, dt):
    scaling = config.scalingFactor
    grad = np.zeros(temp.shape)

    offDiag = dt*(N*N)*scaling
    diag = - 4*dt*(N*N)*scaling

    for row in range(1, 2*N+2):
        for col in range(1, N+1):
            grad[row][col] = diag*temp[row][col] + offDiag*temp[row+1][col] + offDiag*temp[row-1][col] + offDiag*temp[row][col+1] + offDiag*temp[row][col-1]

    temp = temp + grad
    return temp