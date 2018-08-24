import numpy as np
import matplotlib.pyplot as plt


betas = [0.0001, 0.001, 0.01, .1, 1, 10, 100, 1000]
angles = np.linspace(-np.pi/2, np.pi/2, 1000)
#angles = np.linspace(0, 2, 1000)
plt.figure()
for beta in betas:
    phi_s = 5.275 * beta**(-0.3) * angles**2

    plt.plot(angles, phi_s, label=str(beta))


plt.xlabel('slope angle (radians)')
plt.ylabel('$\Phi_s$')
plt.legend(title='beta value')
plt.ylim(-5, 5)
plt.show()
