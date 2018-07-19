import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

'''
0 = Cell is unburnable [U].
1 = Cell is flammable, but not currently ignited [N].
2 = Cell is flammable and is ignited, but fuel is not yet consumed [I].
3 = All fuel in cell has been consumed by the fire [C].

'''

def init_grid(len_side):
    grid = np.ones((len_side, len_side))

    return grid


def rate_of_spread_xy(R0, Ebar, alpha, theta, ti):

    '''
    # old formulation:
    ros_x = (1) * R0*((1-Ebar)
                         / (1-Ebar*np.cos(theta-alpha))) * np.cos(-theta)
    ros_y = (1) * R0*((1-Ebar)
                         / (1-Ebar*np.cos(theta-alpha))) * np.sin(-theta)

    return np.abs(ros_x), np.abs(ros_y)
    '''
    # new formulation:
    ros = R0 * (1-Ebar) / (1-Ebar*np.cos(alpha-theta))

    ros_x = ros * np.cos(theta)
    ros_y = ros * np.sin(theta)

    return np.abs(ros_x), np.abs(ros_y)


def rate_of_spread(R0, Ebar, alpha, theta, ti):
    ros = R0 * (1-Ebar) / (1-Ebar*np.cos(alpha-theta))

    return ros


if __name__ == '__main__':

    # ###### cmap #######
    fire_cmap = colors.ListedColormap(['blue', 'limegreen',
                                       'red', 'k'])
    # bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    bounds = [0, 1, 2, 3, 4]
    norm = colors.BoundaryNorm(bounds, fire_cmap.N)
    # ###### cmap #######

    # params
    Ebar = 0.8
    R0 = 0.5

    timesteps = 70
    len_side = 120
    alpha = np.pi

    y_start = 80

    thetas = [0, np.pi/4, np.pi/2, 3*np.pi/4,
              np.pi, 5*np.pi/4, 3*np.pi/2,
              7*np.pi/4]
    bucket_grid = np.zeros((len_side, len_side))
    neighbor_angles = np.arange(0, 2*np.pi, np.pi/4)
    neighbor_locations = np.array([(+1, 0), (+1, +1), (0, +1), (-1, +1),
                                   (-1, 0), (-1, -1), (0, -1), (+1, -1)])
    neighbor_distances = np.array([1, np.sqrt(2),
                                   1, np.sqrt(2),
                                   1, np.sqrt(2),
                                   1, np.sqrt(2)])
    bucket_grid[y_start, len_side/2] = 1.

    Ebar_grid = Ebar*np.ones((len_side, len_side))
    R0_grid = R0*np.ones((len_side, len_side))
    # Ebar_grid = np.random.random((len_side, len_side))
    # R0_grid = np.random.random((len_side, len_side))

    state_grid = np.ones_like(bucket_grid)
    state_grid[y_start, len_side/2] = 2

    plt.figure()
    plt.imshow(state_grid, cmap=fire_cmap, norm=norm)
    plt.axis('off')
    plt.colorbar()
    plt.savefig('../plots/im000.png')

    for ti in range(timesteps):
        print 'working on timestep', ti
        I_inds_x, I_inds_y = np.where(state_grid == 2)
        for n in range(len(I_inds_x)):
            if np.all(state_grid[I_inds_x[n]-1:I_inds_x[n]+2,
                                 I_inds_y[n]-1:I_inds_y[n]+2] != 1):
                state_grid[I_inds_x[n], I_inds_y[n]] = 3

            else:
                ''' OLD WAY
                for neighbor_num, theta in enumerate(neighbor_angles):
                    ros = rate_of_spread(R0_grid[I_inds_x[n], I_inds_y[n]],
                                         Ebar_grid[I_inds_x[n], I_inds_y[n]],
                                         alpha, theta, ti)
                    bucket_grid[I_inds_x[n]
                                + neighbor_locations[neighbor_num, 0],
                                I_inds_y[n]
                                + neighbor_locations[neighbor_num, 1]] \
                        += (ros / neighbor_distances[neighbor_num])
                    new_I_inds_x, new_I_inds_y = np.where((bucket_grid >= 1.)
                                                          & (state_grid != 3))
                    state_grid[new_I_inds_x, new_I_inds_y] = 2
                '''

                # NEW WAY, vectorized
                neighbor_nums = np.arange(0, len(neighbor_angles))
                ros = rate_of_spread(R0_grid[I_inds_x[n], I_inds_y[n]],
                                     Ebar_grid[I_inds_x[n], I_inds_y[n]],
                                     alpha, neighbor_angles, ti)
                bucket_grid[I_inds_x[n]
                            + neighbor_locations[neighbor_nums, 0],
                            I_inds_y[n]
                            + neighbor_locations[neighbor_nums, 1]] \
                    += (ros / neighbor_distances[neighbor_nums])
                new_I_inds_x, new_I_inds_y = np.where((bucket_grid >= 1.)
                                                      & (state_grid != 3))
                state_grid[new_I_inds_x, new_I_inds_y] = 2

        plt.figure()
        plt.imshow(state_grid, cmap=fire_cmap, norm=norm)
        plt.axis('off')
        plt.colorbar()
        plt.savefig('../plots/im%s.png' % str(ti+1).zfill(3))
        plt.close()
