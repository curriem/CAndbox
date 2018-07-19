import numpy as np
import matplotlib.pyplot as plt


def ros(theta, alpha):
    R0 = 1
    ecc = 0.7
    R = R0 * (1-ecc) / (1-ecc*np.cos(alpha - theta))
    return R


ax = plt.subplot(111, projection='polar')

thetas = np.linspace(0, 2*np.pi, 100)
for a in np.arange(0, 2*np.pi, np.pi/4):
    Rs = ros(thetas, a)

    ax.plot(thetas, Rs, label=str(a))
    for c in np.arange(0, 2*np.pi, np.pi/4):
        ax.scatter(c, ros(c, a))

ax.axhline(0)
ax.axvline(0)
plt.show()
