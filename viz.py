import numpy as np
import matplotlib.pyplot as plt

e = np.loadtxt("fitness.dat")
plt.plot(e)
plt.xlabel("Generations")
plt.ylabel("Fitness")
plt.show()
