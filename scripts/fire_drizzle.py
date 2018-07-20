import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

'''
0 = Cell is unburnable [U].
1 = Cell is flammable, but not currently ignited [N].
2 = Cell is flammable and is ignited, but fuel is not yet consumed [I].
3 = All fuel in cell has been consumed by the fire [C].

'''

correction = True

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


def rate_of_spread(R0, U_eq, alpha, phi, ti):

    LW = 0.936*np.exp(50.5*U_eq) + 0.461*np.exp(-30.5*U_eq) - 0.397
    Ebar = np.sqrt(1 - 1/LW**2)
    theta = phi - alpha

    if correction:

        # corrected rate of spread
        k1 = 1
        k2 = 2
        k3 = 0
        k4 = 0
        k5 = 0

        R0_prime = k1*R0
        Ebar_prime = np.sqrt(1 - 1/LW**k2)
        if (k3 == 0) and (k4 == 0):
            theta_prime = theta
        else:
            if np.abs(theta) <= np.pi/2:
                theta_prime = np.arctan(np.tan(theta) / LW**k3)
            elif np.abs(theta) > np.pi/2:
                theta_prime = np.arctan(np.tan(theta) / LW**k4)

            else:
                assert False, 'weird theta'
        ros = (R0_prime * (1-Ebar_prime) / (1-Ebar_prime*np.cos(theta_prime))
               - R0_prime*((1-Ebar_prime)/(1+Ebar_prime)
                           - (1-Ebar)/(1+Ebar))
               * (np.abs(theta_prime)/np.pi)**k5)

        assert Ebar_prime == Ebar, '%.03f =/= %.03f' % (Ebar_prime, Ebar)
        assert R0_prime == R0, '%.03f =/= %.03f' % (R0_prime, R0)
        assert theta_prime == theta, '%.03f =/= %.03f' % (theta_prime, theta)

    else:
        ros = R0 * (1-Ebar) / (1-Ebar*np.cos(theta))

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
    U_eq = 0.012
    R0 = 0.5

    timesteps = 70
    len_side = 120
    alpha = np.pi - np.pi/4

    y_start = 80

    bucket_grid = np.zeros((len_side, len_side))
    neighbor_angles = np.arange(0, 2*np.pi, np.pi/4)
    neighbor_locations = np.array([(+1, 0), (+1, +1), (0, +1), (-1, +1),
                                   (-1, 0), (-1, -1), (0, -1), (+1, -1)])
    neighbor_distances = np.array([1, np.sqrt(2),
                                   1, np.sqrt(2),
                                   1, np.sqrt(2),
                                   1, np.sqrt(2)])
    bucket_grid[y_start, len_side/2] = 1.

    Ueq_grid = U_eq*np.ones((len_side, len_side))
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
                # OLD WAY
                for neighbor_num, phi in enumerate(neighbor_angles):
                    ros = rate_of_spread(R0_grid[I_inds_x[n], I_inds_y[n]],
                                         Ueq_grid[I_inds_x[n], I_inds_y[n]],
                                         alpha, phi, ti)
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
                                     Ueq_grid[I_inds_x[n], I_inds_y[n]],
                                     alpha, neighbor_angles, ti)
                bucket_grid[I_inds_x[n]
                            + neighbor_locations[neighbor_nums, 0],
                            I_inds_y[n]
                            + neighbor_locations[neighbor_nums, 1]] \
                    += (ros / neighbor_distances[neighbor_nums])
                new_I_inds_x, new_I_inds_y = np.where((bucket_grid >= 1.)
                                                      & (state_grid != 3))
                state_grid[new_I_inds_x, new_I_inds_y] = 2
                '''

        plt.figure()
        plt.imshow(state_grid, cmap=fire_cmap, norm=norm)
        plt.axis('off')
        plt.colorbar()
        plt.savefig('../plots/im%s.png' % str(ti+1).zfill(3))
        plt.close()
