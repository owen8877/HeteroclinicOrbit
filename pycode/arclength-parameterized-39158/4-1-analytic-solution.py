import math
import numpy as np
import scipy as sp
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

T = np.linspace(0, -10, 1000)

p0 = (4, 1)
B = np.matrix([[-1, 4], [-4, -1]])
p = np.array([np.dot(sp.linalg.expm(-np.transpose(B) * t), p0) for t in T])
_B = -np.linalg.inv(B + np.transpose(B))
x = np.array([np.dot(_B, np.dot(sp.linalg.expm(-np.transpose(B) * t), p0)).A1 for t in T])

plt.figure(1)
plt.plot(p[:, 0], p[:, 1])
plt.figure(2)
plt.plot(x[:, 0], x[:, 1])
plt.show()