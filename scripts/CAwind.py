#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 13:00:35 2017

@author: mcurrie


key:
-1 = has been on fire and can no longer ignite
0 = not on fire and can ignite at a later time
1 = on fire for one time step
2 = on fire for two time steps
3 = on fire for three time steps
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_frame(frame, vmin=-1, vmax=1, title=None):
    plt.figure()
    plt.imshow(frame, interpolation='nearest',
               vmin=vmin, vmax=vmax, cmap='jet')
    plt.axis('off')
    plt.colorbar()
    plt.title(title)
    plt.close()


def probMatrix(numLayers, masterProb, windMag, windDir):
    matDim = numLayers*2 + 1  # the dimensions of the matrix of neighbors

    # Compute the distance scales
    mat = np.empty((matDim, matDim))
    center = [numLayers, numLayers]
    for i in range(matDim):
        for j in range(matDim):
            mat[i, j] = np.linalg.norm(np.array([i, j])
                                       - np.array(center))
    mat[center[0], center[1]] = np.inf
    mat = 1/mat
    mat *= masterProb

    # Add wind
    windMat = np.empty_like(mat)
    for i in range(matDim):
        for j in range(matDim):
            if i == j == numLayers:
                windMat[i, j] = 0
                continue
            pixVector = np.array([-(i - center[0]), (j - center[1])])
            ThetaC = (np.arctan2((pixVector[1]), (pixVector[0]+0.000001)))
            cosAlpha = np.cos((windDir - ThetaC))
            windMat[i, j] = cosAlpha
    windMat *= windMag

    probMat = mat + windMat
    probMat[center[0], center[1]] = 0

    negInds = np.array(np.where(probMat < 0))
    probMat[negInds[0], negInds[1]] = 0

    return probMat


def randProbMatrix(numLayers):
    matDim = 2*numLayers + 1
    randMatrix = np.random.rand(matDim, matDim)
    return randMatrix


def burnTimeGrid(N, master, homogeneous=False):
    grid = np.ones((N, N))
    grid[N/3:2*N/3, :] = 2
    grid[2*N/3:, :] = 3
    if homogeneous:
        grid = master*np.ones((N, N))

    return grid


def initialFire(N, fireType):

    fireGrid = np.zeros((N, N))

    if fireType == 'line':
        fireGrid[9*N/10, N/3:2*N/3] = 1
        # fireGrid[N/2+1,2*N/5:3*N/5]=1
    elif fireType == 'point':
        fireGrid[N/2, N/2-1:N/2+2] = 1
        fireGrid[N/2-1, N/2-1:N/2+2] = 1
        fireGrid[N/2+1, N/2-1:N/2+2] = 1
    elif fireType == 'parabola':
        x = np.arange(N/3, 2*N/3)
        y = 0.01*(x-N/2+1)**2 + N/2
        fireGrid[np.rint(y).astype('int'),
                 np.rint(x).astype('int')] = 1
    else:
        print 'Specified fire type not recognized, doing point.'
        fireGrid[N/2, N/2] = 1

    return fireGrid


def main(bigMati, leng, tstep):

    plotting = False   # plot at the end?
    save = False       # save the plots?
    N = leng   # length of a side
    eps = 10**(-3)
    fireType = 'line'
    numLayers = 1
    masterBurnTime = 3
    numTimesteps = tstep

    savePath = '/home/dr17c/FireDynamics/data/'

    fireGrid = initialFire(N, 'line')

    fireTimeseries = [np.copy(fireGrid)]    # initiate a timeseries for later

    masterProb = 0.15*0.75         # Adjusts Pmaster
    alpha = 1.             # Sink strength
    beta = 5.               # Background Flow Strength
    gamma = 0.75               # Tanh Coefficient
    uD = np.zeros((N, N))
    vD = np.zeros((N, N))

    # Uniform flow (Upwards)
    def phi(x, y):
        return (1.*y)

    burnTimes = burnTimeGrid(N, masterBurnTime, homogeneous=True)
    for t in range(numTimesteps):
        print 'working on timestep %i' % t
        fireInds = np.array(np.where(fireGrid > 0))
        for [i, j] in fireInds.T:
            uD[i, j] = beta
            vD[i, j] = 0.

        for [x, y] in fireInds.T:
            centerD = np.array([x, y])

            for [i, j] in fireInds.T:
                pntD = np.array([i, j])
                rx = pntD[0] - centerD[0]
                ry = pntD[1] - centerD[1]
                r2 = rx**2 + ry**2
                uD[i, j] = uD[i, j] - (alpha/(2*np.pi))*rx/(r2+eps**2.)
                vD[i, j] = vD[i, j] - (alpha/(2*np.pi))*ry/(r2+eps**2.)

        for [x, y] in fireInds.T:

            neighbors = fireGrid[x-numLayers:x+numLayers+1,
                                 y-numLayers:y+numLayers+1]
            if np.array_equal(np.array(neighbors.shape),
                              np.array([2.*numLayers+1,
                                        2.*numLayers+1])) == False:
                continue
            randProb = randProbMatrix(numLayers)

            AlphaD = np.arctan2(vD[x, y],
                                (uD[x, y]+eps))

            windMag = np.tanh(gamma*np.linalg.norm([uD[x, y], vD[x, y]]))

            windMag *= masterProb

            probThresh = probMatrix(numLayers, masterProb, windMag, AlphaD)

            probThresh *= (np.ones_like(probThresh)
                           - np.abs(neighbors))

            newFireInds = np.array(np.where(randProb < probThresh))

            neighbors[newFireInds[0], newFireInds[1]] = 1

            fireGrid[x-numLayers:x+numLayers+1,
                     y-numLayers:y+numLayers+1] = neighbors
            fireGrid[x, y] += 1

        burntInds = np.array(np.where(fireGrid > burnTimes))
        fireGrid[burntInds[0], burntInds[1]] = -1
        plot_frame(fireGrid, vmin=-1, vmax=4, title=t+1)
        burningInds = np.array(np.where(fireGrid > 0))
        bigMat[burningInds[0], burningInds[1], t] += 1
        # have to do this because python is stupid sometimes
        temp = np.copy(fireGrid)
        fireTimeseries.append(temp)

    fireTimeseries = np.array(fireTimeseries)

    if save:
        np.save(savePath+'CA_type=%s_layers=%i_burnTime=%i_wind=%s.npy'
                % (fireType, numLayers,
                   masterBurnTime, str(windMag)), fireTimeseries)

    if plotting:
        for frame in fireTimeseries:
            plt.figure()
            plt.imshow(frame, interpolation='nearest',
                       cmap='jet', vmin=-1, vmax=1)
            plt.axis('off')
            plt.colorbar()

    return fireGrid

# main()

#"""
numreps = 1
tstep = 70  # Number of time steps
N = 100 # Grid Size
bigMat = np.zeros([N,N,tstep])
for i in range(numreps):
    main(bigMat,N,tstep)

np.save('frames',bigMat)
#"""
for i in range(tstep):
    frame = bigMat[:,:,i]
    plt.figure(i)
    plt.imshow(frame, interpolation='nearest', cmap='jet', vmin=2, vmax=3)
    plt.axis('off')
    plt.colorbar()
    plt.title(['Timestep:', i+1])
