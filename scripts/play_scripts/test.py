import numpy as np
import matplotlib.pyplot as plt

x, y = np.meshgrid(np.linspace(-100,100,200), np.linspace(-100,100,200))
d = np.sqrt(x*x+y*y)
sigma, mu = 30.0, 1

g = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
g*=10
g = g.astype(int)
print g.shape
print g
#g += 100*np.ones((len_side, len_side), dtype=int)

plt.figure()
plt.imshow(g)
plt.colorbar()
plt.show()
