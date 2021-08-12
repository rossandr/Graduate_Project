import numpy as np
from scipy.optimize import root
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


####### Know Variables  #######
F_pp = 1/10
emiss = 0.85
sigma = 5.67E-8

####### Temperature of particle assumption ########
T_p1 = 300 + 273
T_p3 = 25 + 273

#### effective thermal conductivity ######

D = 0.5 * 0.0254
d = np.array([(1/50), (1/30), (1/20), (1/12), (1/8)])
d = d * 0.0254
D_diff = d / D
E1 = 0.3494 + 0.4381 * D_diff  # sato et al. [16]

Temp_c = 300
Temp = 300 + 273

lam_a = 0.026
lam_m = 0.142
lam_w = 0.6
lam_f = 0.554 + 2.24E-3*Temp_c - 9.87E-6*Temp_c**2
lam_da = 0.024 + 7.73E-5*Temp_c - 9.87E-8*Temp_c**2 # for temperatures above 600 deg C

x_w = E1*0.5
x_a = E1*0.95
x_m = 1 - E1

g_a = 1.7
g_b = g_a
g_c = 1-(2*g_a)

k_a = (1/3)*((2/(1+(lam_a/lam_f - 1)*g_a)) + (1/(1 + (lam_a/lam_f -1)*g_c)))
k_w = (1/3)*((2/(1+(lam_w/lam_f - 1)*g_a)) + (1/(1 + (lam_w/lam_f -1)*g_c)))
k_m = (1/3)*((2/(1+(lam_m/lam_f - 1)*g_a)) + (1/(1 + (lam_m/lam_f -1)*g_c)))

K_p = (k_w*x_w*lam_w + k_a*x_a*lam_a + k_m*x_m*lam_m) / (k_w*x_w + k_a*x_a + k_m*x_m)

###### Convective Heat Transfer Coefficient #####
v = 20.93E-6 #m^2/s
alpha = 29.9E-6 #m^2/s
Ts = 300 #Kelvin
T_inf = 25 + 273 #Kelvin
g = 9.81 #m/s^2
Beta = 3.12E-3 #1/K
R = np.array([(1/50), (1/30), (1/20), (1/12), (1/8)])
R = R*0.0254
Pr = 0.700
K_air = 30.0E-3 #(W/m*K)


Ra_D = (g*Beta*(Ts - T_inf)*R**3)/(v*alpha)
Nu_d = 2 + ((0.586*Ra_D**(1/4))/((1 + (0.469/Pr)**(9/16))**(4/9)))

h_conv = (Nu_d * K_air)/R


#### Heat transfer portions #####

# Rad_in = F_pp * emiss * sigma * (T_p1**4 - T_p2**4)
# Rad_out = F_pp * emiss * sigma * (T_p2**4 - T_p3**4)
# Cond_in = (K_p/L)*(T_p1 - T_p2)
# Cond_out = (K_p/L)*(T_p2 - T_p3)
# Conv_out = h_conv * (T_p2 - T_f)


######### Solve for equation to be zero ##########
def solve(x, K_air, h_conv, K_p):

    R = np.array([(1 / 50), (1 / 30), (1 / 20), (1 / 12), (1 / 8)])
    R = R * 0.0254
    # print(R[0])
    R1 = R[4]
    h_conv = h_conv[4]
    K_p = K_p[4]

    # print(h_conv)
    # print(K_p)

    t = x[0]
    # t = np.arange(0.0000001, 3, 0.00000005)
    # # t = 1
    Tf = 25 + 273
    Tr = 300 + 273
    kf = K_air
    cp_p = 1.76
    rho_p = 112  # kg/m^3
    Tp_1 = 300 + 273
    Tp = 250 + 273

    return (-Tp + ((F_pp*emiss*sigma*Tr**4 + (K_p/R1)*Tp_1 + h_conv*Tf - np.exp(((-3*t)/(R1*rho_p*cp_p)*((K_p/R1)+h_conv))))/((K_p/R1)+h_conv)))

x = [0.00001]

sol = root(solve, x, args=(K_air, h_conv, K_p,))

print('time:', sol.x[0])


########### Test ###############

# R = np.array([(1 / 50), (1 / 30), (1 / 20), (1 / 12), (1 / 8)])
# R = R * 0.0254
# # print(R[0])
# R1 = R[2]
# h_conv = h_conv[2]
# K_p = K_p[2]
#
# # print(h_conv)
# # print(K_p)
#
# t = np.arange(0.0000001, 3, 0.00000005)
# # t = 1
# Tf = 25 + 273
# Tr = 300 + 273
# kf = K_air
# cp_p = 1.76
# rho_p = 112 #kg/m^3
# Tp_1 = 350+273
# # Tp = 250 + 273
#
# # Tp = Tf + ((emiss*sigma*Tr**4*R)/(kf))*(1 - np.exp((-3*kf*t)/(rho_p*cp_p*R**2)))
# Tp = (F_pp*emiss*sigma*Tr**4 + (K_p/R1)*Tp_1 + h_conv*Tf - np.exp(((-3*t)/(R1*rho_p*cp_p)*((K_p/R1)+h_conv))))/((K_p/R1)+h_conv)
#
# print(Tp - 273)
#
# plt.plot(t, Tp)
# plt.show()


