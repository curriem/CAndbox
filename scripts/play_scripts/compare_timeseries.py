import numpy as np
from matplotlib import colors

import matplotlib.pyplot as plt
# ###### cmap #######
fire_cmap = colors.ListedColormap(['blue', 'limegreen',
                                   'red', 'k'])
bounds = [0, 1, 2, 3, 4]
norm = colors.BoundaryNorm(bounds, fire_cmap.N)
# ###### cmap #######

ts1 = np.load('../model_timeseries/mountain.npy')
ts2 = np.load('../model_timeseries/nomountain.npy')

print ts1[100]
n, x, y = ts1.shape

for m in range(n):
    plt.figure()
    im = np.hstack((ts1[m], ts2[m]))
    plt.axis('off')
    plt.imshow(im, cmap=fire_cmap, norm=norm)
    plt.savefig('../plots/compare%s.png' % str(m).zfill(3))
    plt.close()
