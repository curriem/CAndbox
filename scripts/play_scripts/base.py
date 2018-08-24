import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


'''
0 = Cell is unburnable [U].
1 = Cell is flammable, but not currently ignited [N].
2 = Cell is flammable and is ignited, but fuel is not yet consumed [I].
3 = All fuel in cell has been consumed by the fire [C].

'''


def init_grid(N_side):

    grid = np.ones((N_side, N_side))

    return grid


def rate_of_spread(ecc, R0, theta):
    # rate of spread in cells/timestep
    ros = R0 * (1 - ecc)/(1 - ecc*np.cos(theta))

    return ros


def neighbor_angle(max_spread_angle, x, y):

    '''
    -------------------------
    |       |       |       |
    | (2,0) | (2,1) | (2,2) |
    |       |       |       |
    -------------------------
    |       |       |       |
    | (1,0) | (1,1) | (1,2) |
    |       |       |       |
    -------------------------
    |       |       |       |
    | (0,0) | (0,1) | (0,2) |
    |       |       |       |
    -------------------------
    '''
    # layers: angle, distance
    neighbors = np.empty((3, 3))

    neighbors[1, 2] = 0 - max_spread_angle
    neighbors[2, 2] = np.pi/4 - max_spread_angle
    neighbors[2, 1] = np.pi/2 - max_spread_angle
    neighbors[2, 0] = 3*np.pi/4 - max_spread_angle
    neighbors[1, 0] = np.pi - max_spread_angle
    neighbors[0, 0] = 5*np.pi/4 - max_spread_angle
    neighbors[0, 1] = 3*np.pi/2 - max_spread_angle
    neighbors[0, 2] = 7*np.pi/4 - max_spread_angle

    neighbor_angle = neighbors[x, y]
    return neighbor_angle


if __name__ == '__main__':

    # params
    ecc = 0.7   # eccentricity of the fire ellipse
    R0 = 1      # rate of spread with no wind or slope
    # theta is the direction of spread measured from the angle of max spread
    timesteps = 50
    N_side = 200
    max_spread_angle = np.pi/2

    # ###### cmap #######
    fire_cmap = colors.ListedColormap(['blue', 'limegreen',
                                       'red', 'k'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = colors.BoundaryNorm(bounds, fire_cmap.N)
    # ###### cmap #######

    grid = np.empty((0, N_side, N_side))
    initial_grid = init_grid(N_side)
    initial_grid[N_side/2, N_side/2] = 2
    initial_grid[2*N_side/3 - 2:2*N_side/3 + 3,
                 N_side/2 - 2:N_side/2 + 3] = 0
    initial_grid[2*(N_side)/3 - 2:2*(N_side)/3 + 3,
                 (N_side-30)/2 - 2:(N_side-30)/2 + 3] = 0
    initial_grid[2*(N_side)/3 - 2:2*(N_side)/3 + 3,
                 (N_side+30)/2 - 2:(N_side+30)/2 + 3] = 0
    grid = np.concatenate((grid, [initial_grid]))
    for t in range(timesteps):
        old_step = grid[-1, :, :]
        new_step = old_step.copy()
        I_inds_x, I_inds_y = np.where(old_step == 2)
        print 'Timestep %i: %i pixels on fire' % (t, len(I_inds_x))
        for n in range(len(I_inds_x)):
            neighbor_grid = old_step[I_inds_x[n]-1:I_inds_x[n]+2,
                                     I_inds_y[n]-1:I_inds_y[n]+2]
            if np.all(neighbor_grid != 1):
                # if all cells in neighborhood are I or C, set cell to C
                new_step[I_inds_x[n],
                         I_inds_y[n]] = 3
            else:
                # if a cell is labeled I and one or more neighbors are N
                N_inds_x, N_inds_y = np.where(neighbor_grid == 1)
                for m in range(len(N_inds_x)):
                    # angle wrt max fire spread angle
                    N_angle = neighbor_angle(max_spread_angle,
                                             N_inds_x[m],
                                             N_inds_y[m])
                    ros = rate_of_spread(ecc, R0, N_angle)
                    if round(ros) == 1.:
                        new_step[I_inds_x[n] + N_inds_x[m] - 1,
                                 I_inds_y[n] + N_inds_y[m] - 1] = 2
                    else:
                        pass
        grid = np.concatenate((grid, [new_step]))
        ts, i, j = grid.shape

        plt.figure()
        plt.imshow(grid[t, :, :], cmap=fire_cmap, norm=norm,
                   origin='lower')
        plt.axis('off')
        plt.savefig('../plots/step%s.png' % str(t).zfill(3))
        plt.close()
