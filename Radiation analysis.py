import numpy as np
from matplotlib import pyplot as plt

# Particle Properties
K_s = 0.142
Cp_s = 1733  # k/kg*k
rho_s = 660  # kg/m^3
# Fluid Properties
K_f = 26.3E-3
d = np.array([(1 / 8), (1 / 12), (1 / 20), (1 / 30), (1 / 50)])
d_m = d * 0.0254
# d = np.array([(1/50), (1/30), (1/20), (1/12), (1/8)])
# d = d1[::-1]
D = 0.1
D_accurate = D/d_m
D_diff = d / D
# porosity
E1 = 0.3494 + 0.4381 * D_diff  # sato et al. [16]
print(E1)
E2 = 0.681 - 1.363 * D_diff + 2.241 * (D_diff ** 2)
E3 = 0.4 + 0.05 * D_diff + 0.412 * (D_diff ** 2)
E = np.array([E1, E2, E3])
# Thermal Conductivity
Keff_min = 1 / (((1 - E) / K_s) + E / K_f)
Keff_max = E * K_f + (1 - E) * K_s
K_small_p = ((K_f + 2 * K_s - 2 * E * (K_s - K_f)) / (K_f + 2 * K_s + E * (K_s - K_f))) * K_s
K = np.array([Keff_min, Keff_max, K_small_p])
# print(Keff_min)
# print(Keff_max)
# print(K_small_p)
### Heat Rate
L = 0.2
h = 0.1
A = L * h
T1 = 300
T2 = 20

fig1, ax1 = plt.subplots()

plt.plot(d, Keff_min[1], label='Keff_min')
plt.plot(d, Keff_max[1], label='Keff_max')
plt.plot(d, K_small_p[1], label='K for small particles')
ax1.set_xlabel('Particle Diameter [in]')
ax1.set_ylabel('Thermal Conductivity [W/m*k]')
ax1.set_title('Thermal Conductivity vs. Particle Size')
ax1.legend()
ax1.set_xlim([0.125, 0.02])
# plt.show()
print(K_small_p)

q = ((K * A) / L) * (T1 - T2)
fig2, ax2 = plt.subplots()

plt.plot(d, q[0][1], label='Keff_min')
plt.plot(d, q[1][1], label='Keff_max')
plt.plot(d, q[2][1], label='K for small particles')
ax2.set_xlabel('Particle Diameter [in]')
ax2.set_ylabel('Heat Rate [W]')
ax2.set_title('Heat Rate vs. Particle Size')
ax2.set_xlim([0.125, 0.02])
ax2.legend()
plt.show()
########### Convective Heat Transfer Coefficient ############
# A_b = 0.2032**2
A_b = np.pi*(0.127/2)**2
rho_f = 0.8711
mu_f = 230.1E-7
Pr_f = 0.690
cp_f = 1.014E+3
K_f = 33.8E-3
sccm = 100
Q = sccm*0.000000017
m_p_s = Q/A_b
print(d)
Re_p = (d_m*m_p_s*rho_f)/(mu_f*(1-E[0]))
print(Re_p)
# colburn factors
Ar_m = ((d_m**3)*9.81*rho_f*(rho_s - rho_f)*(1 - E[0])**2)/mu_f**2
print(Ar_m)
E_ih3 = 0.018*(Ar_m/Re_p**2)**0.25  # Bhattacharyya
E_ih1 = 2.52 * Re_p**(-0.5)
E_ih2 = 1.25*Re_p**0.58
Ejh = np.array([E_ih1, E_ih2, E_ih3])
Nu = Ejh * Re_p * Pr_f**(2/3)
h_jh = (Nu*K_f)/d_m
print(E_ih1)
h1 = (2.06*rho_f*m_p_s*cp_f)/(E[0]*Pr_f**(2/3)*Re_p**0.575)
Nu_kta = 1.27*(Pr_f**(1/2)/E[0]**1.28)*Re_p**0.36 + 0.033*(Pr_f**(1/2)/E[0]**1.07)*Pr_f**0.86
h2 = (Nu_kta*K_f)/d_m
print(h2)

# fig3, ax3 = plt.subplots()
#
# ax3.plot(d, h1, label='Correlation 1')
# ax3.plot(d, h2, label='Correlation 2')
# # plt.plot(d, h_jh[2], label='Bhattacharyya')
# ax3.set_xlabel('Particle Diameter [in]')
# ax3.set_ylabel('Convective Heat Transfer Coefficient [W/m^2*K]')
# ax3.set_title('Convective Heat Transfer Coefficient vs. Particle Size')
# ax3.set_xlim([0.125, 0.02])
# ax3.legend()
#
#
# A_sphere = 4*np.pi*(d_m/2)**2
# q_jh = h_jh*A_sphere*(250-20)
#
# q1 = h1*A_sphere*(250 - 20)
# q2 = h2*A_sphere*(250 - 20)
# #
# # plt.figure(4)
#
# fig4, ax4 = plt.subplots()
#
# ax4.plot(d, q1, label='Correlation 1')
# ax4.plot(d, q2, label='Correlation 2')
# # ax4.plot(d, q_jh[2], label='Bhattacharyya')
#
# ax4.set_xlabel('Particle Diameter [in]')
# ax4.set_ylabel('q Heat Rate [W]')
# ax4.set_title('Convective Heat Transfer Coefficient vs. Particle Size')
# ax4.legend()
#
# ax4.set_xlim([0.125, 0.02])
# print()
# print('colburn relationships')
# print(h_jh)
# plt.show()



