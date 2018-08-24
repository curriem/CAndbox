import numpy as np
import matplotlib.pyplot as plt


U_eq = np.arange(0, 10, 0.01)
LW = 0.936*np.exp(50.5*U_eq) + 0.461*np.exp(-30.5*U_eq) - 0.397
U_eq = 0.78
LW = 0.936*np.exp(0.2566*U_eq) + 0.461*np.exp(-0.1548*U_eq) - 0.397

Ebar = np.sqrt(1 - 1/LW**2)
print Ebar
assert False
ind = np.where(np.isclose(Ebar, 0.5, atol=0.01))
print Ebar[ind], LW[ind]
plt.figure()
plt.plot(U_eq, Ebar)
plt.xlabel('U_eq')
plt.ylabel('Ebar')
plt.show()
