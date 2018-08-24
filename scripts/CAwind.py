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
import argparse

parser = argparse.ArgumentParser(description='CA wind')
parser.add_argument('save_file',
                    help='name of file you want to save output, (.npy)')
args = parser.parse_args()


def plot_frame(frame, vmin=-1, vmax=1, title=None):
    plt.figure()
    plt.imshow(frame, interpolation='nearest',
               vmin=vmin, vmax=vmax, cmap='jet')
    plt.axis('off')
    plt.colorbar()
    plt.title(title)
    plt.close()


def probability_matrix(neighbor_radius, masterProb, windMag, windDir):
    matDim = neighbor_radius*2 + 1  # the dimensions of the matrix of neighbors

    # Compute the distance scales
    mat = np.empty((matDim, matDim))
    center = [neighbor_radius, neighbor_radius]
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
            if i == j == neighbor_radius:
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


def rand_probs(neighbor_radius):
    matDim = 2*neighbor_radius + 1
    randMatrix = np.random.rand(matDim, matDim)
    return randMatrix


def burnTimeGrid(len_side, master, homogeneous=False):
    grid = np.ones((len_side, len_side))
    grid[len_side/3:2*len_side/3, :] = 2
    grid[2*len_side/3:, :] = 3
    if homogeneous:
        grid = master*np.ones((len_side, len_side))

    return grid


def initialFire(len_side, fireType):

    fireGrid = np.zeros((len_side, len_side))

    if fireType == 'line':
        fireGrid[9*len_side/10, len_side/3:2*len_side/3] = 1
        # fireGrid[len_side/2+1,2*len_side/5:3*len_side/5]=1
    elif fireType == 'point':
        fireGrid[len_side/2, len_side/2-1:len_side/2+2] = 1
        fireGrid[len_side/2-1, len_side/2-1:len_side/2+2] = 1
        fireGrid[len_side/2+1, len_side/2-1:len_side/2+2] = 1
    elif fireType == 'parabola':
        x = np.arange(len_side/3, 2*len_side/3)
        y = 0.01*(x-len_side/2+1)**2 + len_side/2
        fireGrid[np.rint(y).astype('int'),
                 np.rint(x).astype('int')] = 1
    else:
        print 'Specified fire type not recognized, doing point.'
        fireGrid[len_side/2, len_side/2] = 1

    return fireGrid


def make_timeseries(frames, len_side, num_t):

    epsilon = 1e-3
    fireType = 'line'
    neighbor_radius = 1
    masterBurnTime = 3

    fireGrid = initialFire(len_side, 'line')

    fireTimeseries = [np.copy(fireGrid)]    # initiate a timeseries for later

    masterProb = 0.15*0.75         # Adjusts Pmaster
    alpha = 1.             # Sink strength
    beta = 5.               # Background Flow Strength
    gamma = 0.75               # Tanh Coefficient
    uD = np.zeros((len_side, len_side))
    vD = np.zeros((len_side, len_side))

    burnTimes = burnTimeGrid(len_side, masterBurnTime, homogeneous=True)
    for t in range(num_t):
        print 'working on timestep %i' % t
        fire_x, fire_y = np.array(np.where(fireGrid > 0))

        uD[fire_x, fire_y] = beta
        vD[fire_x, fire_y] = 0.

        for n in range(len(fire_x)):
            rx = fire_x - fire_x[n]
            ry = fire_y - fire_y[n]
            r2 = rx**2 + ry**2  # is this supposed to be normalized?
            offset_term = (alpha/(2*np.pi))*rx/(r2+epsilon**2.)
            uD[fire_x, fire_y] = uD[fire_x, fire_y] - offset_term
            vD[fire_x, fire_y] = vD[fire_x, fire_y] - offset_term

        for n in range(len(fire_x)):
            neighbors = fireGrid[fire_x[n]-neighbor_radius:
                                 fire_x[n]+neighbor_radius+1,
                                 fire_y[n]-neighbor_radius:
                                 fire_y[n]+neighbor_radius+1]

            try:
                assert neighbors.shape == (2.*neighbor_radius+1,
                                           2.*neighbor_radius+1)
            except AssertionError:
                continue

            randProb = rand_probs(neighbor_radius)

            AlphaD = np.arctan2(vD[fire_x[n], fire_y[n]],
                                (uD[fire_x[n], fire_y[n]]))

            windMag = np.tanh(gamma*np.linalg.norm([uD[fire_x[n],
                                                       fire_y[n]],
                                                    vD[fire_x[n],
                                                       fire_y[n]]]))

            windMag *= masterProb

            probThresh = probability_matrix(neighbor_radius,
                                            masterProb, windMag,
                                            AlphaD)

            probThresh *= (np.ones_like(probThresh)
                           - np.abs(neighbors))

            newFireInds = np.array(np.where(randProb < probThresh))

            neighbors[newFireInds[0], newFireInds[1]] = 1

            fireGrid[fire_x[n]-neighbor_radius:
                     fire_x[n]+neighbor_radius+1,
                     fire_y[n]-neighbor_radius:
                     fire_y[n]+neighbor_radius+1] = neighbors

            fireGrid[fire_x[n], fire_y[n]] += 1

        burntInds = np.array(np.where(fireGrid > burnTimes))
        fireGrid[burntInds[0], burntInds[1]] = -1
        plot_frame(fireGrid, vmin=-1, vmax=4, title=t+1)
        burningInds = np.array(np.where(fireGrid > 0))
        frames[burningInds[0], burningInds[1], t] += 1

        temp = np.copy(fireGrid)
        fireTimeseries.append(temp)

    fireTimeseries = np.array(fireTimeseries)

    np.save('../model_timeseries/CA_type=%s_layers=%i_burnTime=%i_wind=%s.npy'
            % (fireType, neighbor_radius,
               masterBurnTime, str(windMag)), fireTimeseries)

    return fireGrid


if __name__ == '__main__':

    num_reps = 1
    t = 70  # Number of time steps
    len_side = 100  # Grid Size
    frames = np.zeros([len_side, len_side, t])
    for i in range(num_reps):
        make_timeseries(frames, len_side, t)

    np.save(args.save_file, frames)
