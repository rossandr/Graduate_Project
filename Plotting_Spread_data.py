import numpy as np
from matplotlib import pyplot as plt

######################## Plotting Spread Rates Vs. Particle Size #####################

x_axis = (
'3.18-2.12', '3.18-2.12',
'2.12-1.27', '2.12-1.27',
'1.27-0.85', '1.27-0.85', '1.27-0.85', '1.27-0.85', '1.27-0.85',
'0.85-0.51',  '0.85-0.51', '0.85-0.51', '0.85-0.51', '0.85-0.51', '0.85-0.51',
'<0.51', '<0.51', '<0.51')

y_axis = ((0, 0), (0, 0), (1.258, 1.12, 0, 0), (1.38, 1.341, 1.37, 1.38, 1.361, 0.83), (1.6365, 1.58))

x_axis1 = (1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5)
y_axis1 = (0, 0,
           0, 0,
           1.158, 1.12, 0, 0,  0.945,
           1.38, 1.341, 1.37, 1.38, 1.361, 1.326,
           1.6365, 1.58, 1.5018,)

flame_data_x1 =(
'3.18-2.12', '3.18-2.12', '3.18-2.12',
'2.12-1.27', '2.12-1.27',
'1.27-0.85',
'0.85-0.51',
'<0.51', '<0.51')
flame_data_y1 = (0, 0, 0,
                 0, 0,
                 1.09,
                 1.31,
                 1.59, 1.536)

low_temp_data_x1 = ('3.18-2.12', '3.18-2.12',
                    '2.12-1.27',
                    '1.27-0.85', '1.27-0.85',
                    '0.85-0.51',
                    '<0.51', '<0.51')
low_temp_data_y1 = (0, 0,
                    0,
                    0.756, 0.746,
                    1.326,
                    1.426, 1.585)

fig, ax1 = plt.subplots()
ax1.scatter(x_axis[::-1], y_axis1[::-1], marker="x", color='green')
ax1.scatter(flame_data_x1[::-1], flame_data_y1[::-1], marker='o', color='red')
ax1.scatter(low_temp_data_x1[::-1], low_temp_data_y1[::-1], marker='o', color='blue')
ax1.set_xlabel('Particle Size [mm]')
ax1.set_ylabel('Spread rate [mm/min]')
# ax1.set_title('Spread Rate vs. Particle Size')


############# Plotting Spread Rate Vs. Permeability ###############

# y_axis1 = (0, 0,
#            0, 0,
#            1.158, 1.12, 0, 0,  0.945,
#            1.38, 1.341, 1.37, 1.38, 1.361, 1.326,
#            1.6365, 1.58, 1.5018,)

x_axis2 = (0.052, 1.12, 0.521, 0.481, 0.561, 0.723, 0.942, 0.442, 0.412, 0.364, 0.398, 0.551, 0.089, 0.079, 0.084)
y_axis2 = (1.63, 0,  1.258, 1.38, 0.83, 0, 0, 1.341, 1.37, 1.38, 1.361, 1.12, 1.58, 1.5018, 1.326)

flame_data_x2 = (4.55, 4.37, 4.76,
                 1.2, 0.95,
                 0.56,
                 0.43,
                 0.12, 0.112)
flame_data_y2 = (0, 0, 0,
                 0, 0,
                 1.09,
                 1.31,
                 1.59, 1.536)

low_temp_data_x2 = (4.53, 4.27,
                    1.52,
                    0.46, 0.49,
                    0.39,
                    0.14, 0.098)
low_temp_data_y2 = (0, 0,
                    0,
                    0.756, 0.746,
                    1.326,
                    1.426, 1.585)

fig2, ax2 = plt.subplots()
ax2.scatter(x_axis2, y_axis2, marker="x", color='green')
ax2.scatter(flame_data_x2, flame_data_y2, marker='o', color='red')
ax2.scatter(low_temp_data_x2, low_temp_data_y2, marker='o', color='blue')
ax2.set_xlabel('Permeability [um^2]')
ax2.set_ylabel('Spread rate [mm/min]')
# plt.title('Spread Rate vs. Permeability')



################# Plotting Spread Rate Vs. Porosity ####################

x_axis3 = (49.2, 71.1, 65.52, 58.9, 59.08, 65.38, 68, 61.9, 55.4, 58.3, 54.2, 59.3, 58.1)
y_axis3 = (1.64, 0, 1.258, 1.38, 0.83, 0, 0, 1.341, 1.37, 1.38, 1.361, 1.12, 1.58)

flame_data_x3 = (76.2, 71.9, 80.1,
                 71.5, 69.3,
                 63.2,
                 55.6,
                 48.3, 51.2)

flame_data_y3 = (0, 0, 0,
                 0, 0,
                 1.09,
                 1.31,
                 1.59, 1.536)


low_temp_data_x3 = (78.1, 70.1, 75.4,
                    68.6,
                    62.5, 60.3,
                    54.3,
                    48.9, 46.2)

low_temp_data_y3 = (0, 0, 0,
                    0,
                    0.756, 0.746,
                    1.326,
                    1.426, 1.585)

fig3, ax3 = plt.subplots(1, 2, figsize=(10, 6), sharey='row')
ax3[0].scatter(x_axis3, y_axis3, marker="x", color='green')
ax3[0].scatter(flame_data_x3, flame_data_y3, marker='o', color='red')
ax3[0].scatter(low_temp_data_x3, low_temp_data_y3, marker='o', color='blue')
ax3[0].set_xlabel('Porosity [%]')
ax3[0].set_ylabel('Spread rate [mm/min]')
# plt.title('Spread Rate vs. Porosity')


############# Plotting Spread Rate Vs. Permeability ###############

# plt.figure(4)
ax3[1].scatter(x_axis2, y_axis2, marker="x", color='green')
ax3[1].scatter(flame_data_x2, flame_data_y2, marker='o', color='red')
ax3[1].scatter(low_temp_data_x2, low_temp_data_y2, marker='o', color='blue')
ax3[1].set_xlim(0, 1)
ax3[1].set_xlabel('Permeability [um^2]')
# ax3[1].set_ylabel('Spread rate [mm/min]')
# plt.title('Spread Rate vs. Permeability')

for ax in ax3.flat:
    ax.set(ylabel='Spread Rate [mm/min]')

for ax in ax3.flat:
    ax.label_outer()

plt.show()


fig3, (ax3, ax5) = plt.subplots(1, 2, figsize=(10, 6))
ax3.scatter(x_axis3, y_axis3, marker="x", color='green')
ax3.scatter(flame_data_x3, flame_data_y3, marker='o', color='red')
ax3.scatter(low_temp_data_x3, low_temp_data_y3, marker='o', color='blue')
ax3.set_xlabel('Porosity [%]')
ax3.set_ylabel('Spread rate [mm/min]')
# plt.title('Spread Rate vs. Porosity')


############# Plotting Spread Rate Vs. Permeability ###############

# plt.figure(4)
ax5.scatter(x_axis2, y_axis2, marker="x", color='green')
ax5.scatter(flame_data_x2, flame_data_y2, marker='o', color='red')
ax5.scatter(low_temp_data_x2, low_temp_data_y2, marker='o', color='blue')
ax5.set_xlim(0, 1)
ax5.set_xlabel('Permeability [um^2]')
# ax3[1].set_ylabel('Spread rate [mm/min]')
# plt.title('Spread Rate vs. Permeability')
plt.show()


x_axis = (
'3.18-2.12', '3.18-2.12',
'2.12-1.27', '2.12-1.27',
'1.27-0.85', '1.27-0.85', '1.27-0.85', '1.27-0.85', '1.27-0.85',
'0.85-0.51',  '0.85-0.51', '0.85-0.51', '0.85-0.51', '0.85-0.51', '0.85-0.51',
'<0.51', '<0.51', '<0.51')

y_axis1 = ('Non-Smoldering', 'Non-Smoldering',
           'Non-Smoldering', 'Non-Smoldering',
           'Smoldering', 'Smoldering', 'Non-Smoldering', 'Non-Smoldering',  'Smoldering',
           'Smoldering', 'Smoldering', 'Smoldering', 'Smoldering', 'Smoldering', 'Smoldering',
           'Smoldering', 'Smoldering', 'Smoldering',)

fig5, ax5 = plt.subplots()
plt.scatter(x_axis[::-1], y_axis1[::-1], marker="x", color='red')
ax5.set_xlabel('Particle Size [mm]')
ax5.set_ylabel('Smoldering Status')
# ax5.set_title('Spread Rate vs. Particle Size')
ax5.invert_yaxis()
# ax5.set_ylim('Non-Smoldering', 'Smoldering')


plt.show()





