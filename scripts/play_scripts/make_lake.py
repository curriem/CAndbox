import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
grid = np.ones((200, 200))

lake_im = ndimage.imread('polka.png')
lake = np.ones((200, 200))

for n in range(200):
    for m in range(200):
        if np.all(lake_im[n, m, :] == 255):
            pass
        else:
            lake[n, m] = 0

plt.figure()
plt.imshow(lake)
plt.show()
np.save('polka.npy', lake)
