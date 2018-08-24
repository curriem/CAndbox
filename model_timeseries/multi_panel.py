import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt

panel1 = np.load('flat.npy')
panel2 = np.load('wind90.npy')
panel3 = np.load('pyramid.npy')
panel4 = np.load('lake.npy')
print panel1.shape
print panel2.shape
print panel3.shape
print panel4.shape
print np.array_equal(panel1, panel2)
t, x, y = panel1.shape
fire_cmap = colors.ListedColormap(['blue', 'limegreen',
                                   'red', 'k'])
bounds = [0, 1, 2, 3, 4]
norm = colors.BoundaryNorm(bounds, fire_cmap.N)

for n in range(t):
    fig, ax = plt.subplots(2, 2, figsize=(4,4))
    #im = np.hstack((panel1[n, :, :], panel2[n, :, :]))
    #plt.figure()
    ax[0, 0].imshow(panel1[n, :, :], cmap=fire_cmap, norm=norm)
    ax[0, 1].imshow(panel2[n, :, :], cmap=fire_cmap, norm=norm)
    ax[1, 0].imshow(panel3[n, :, :], cmap=fire_cmap, norm=norm)
    ax[1, 1].imshow(panel4[n, :, :], cmap=fire_cmap, norm=norm)
    #plt.imshow(im, cmap=fire_cmap, norm=norm)
    ax[0, 0].axis('off')
    ax[1, 1].axis('off')
    ax[0, 1].axis('off')
    ax[1, 0].axis('off')
    ax[0, 0].set_title('flat terrain')
    ax[0, 1].set_title('N wind to E wind')
    ax[1, 0].set_title('pyramid terrain')
    ax[1, 1].set_title('lake terrain')
    ts = n+1
    ax[0, 0].text(130, -50, 'timestep %i' % ts, 
                 fontsize=18)
    plt.subplots_adjust(wspace=0.25, hspace=0.25)
    plt.savefig('multi_im%s.png' % str(n).zfill(3), bbox_inches='tight')
    plt.close()
