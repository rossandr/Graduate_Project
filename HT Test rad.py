import numpy as np
from scipy.optimize import root
from scipy.optimize import fsolve


R = np.array([(1 / 50), (1 / 30), (1 / 20), (1 / 12), (1 / 8)])
R = R[0] * 0.0254


t = 0.001
Tf = 25 + 273
Tr = 300 + 273
kf = K_air
cp_p = 1.76
rho_p = 112 #kg/m^3
Tp_1 = 300+273
Tp = 250

Tp = Tf + ((emiss*sigma*Tr**4*R)/(kf))*(1 - np.exp((-3*kf*t)/(rho_p*cp_p*R**2)))

print(Tp)

