import math
import numpy as np
import scipy as sp
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

print('Using ArcLength Method')


def xDifferential(x, s, pInt, a, B, C):
    # @param    x - array
    # @return   dx - array
    x = np.transpose(np.matrix(x))
    _p = pInt(s)
    p = np.transpose(np.matrix(_p))

    b = B * x
    bnorm = np.inner(b.A1, (np.linalg.inv(a) * b).A1)

    dx = C / bnorm * (a * p + b)
    return dx.A1

def pDifferential(p, s, xInt, a, B, C):
    # @param    p - array
    # @return   dp - array
    p = np.transpose(np.matrix(p))
    _x = xInt(s)
    x = np.transpose(np.matrix(_x))

    b = B * x
    bnorm = np.inner(b.A1, (np.linalg.inv(a) * b).A1)
    nablab = -B

    dp = -C / bnorm * np.transpose(nablab) * p
    return dp.A1

x0 = (0, 0)
p1 = (4, 1)
sPos = np.linspace(0, 1, 50)
sNeg = np.linspace(1, 0, 50)
xInt = lambda s : np.array([s, 0])
a = np.matrix([[1, 0], [0, 1]])
B = np.matrix([[-1, 4], [-4, -1]])
C = 1

pSol = odeint(pDifferential, p1, sNeg, args=(xInt, a, B, C))
pInt = interp1d(sNeg, np.transpose(pSol))
plt.figure(1)
plt.plot(pSol[:, 0], pSol[:, 1], color='b')

# xSol = odeint(xDifferential, x0, sPos, args=(pInt, a, B, C))
# xInt = interp1d(sPos, np.transpose(xSol))
# plt.figure(2)
# plt.plot(xSol[:, 0], xSol[:, 1], color='b')
plt.show()
