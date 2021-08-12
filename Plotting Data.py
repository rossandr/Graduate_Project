import numpy as np
import matplotlib.pyplot as plt

wind_speed = np.asarray([0, 1, 2, 3, 4, 5])
spread_rate = np.asarray([0.63016, 1.5833, 2.32, 2.596, 3.296, 3.603])
peak_temperature = np.asarray([514.8425, 557.26, 613.35, 634.33, 677.56, 704.83])

spread_uncertainty = np.asarray([0.01, 0.0961, 0.0373, 0.3037, 0.3037, 0.2121])
temperature_uncertainty = np.asarray([23.27, 39.3, 42.38, 83.1, 55.47, 20.79])
windspeed_uncertainty = np.asarray([0.12, 0.12, 0.12, 0.12, 0.12, 0.12])

plt.figure(1, figsize=(7, 7))
plt.errorbar(wind_speed, spread_rate, yerr=spread_uncertainty, xerr=None, capsize=5, fmt='o', capthick=2, ecolor='g')
plt.xlim(0)
plt.ylim([0, 5])
plt.grid()
plt.xlabel('Wind Speed [m/s]')
plt.ylabel('Spread Rate [mm/min]')


plt.figure(2, figsize=(7, 7))
plt.errorbar(wind_speed, peak_temperature, yerr=temperature_uncertainty, xerr=None, capsize=5, fmt='o', capthick=2, ecolor='g')
plt.xlim(0)
# plt.ylim([0, 5])
plt.grid()
plt.xlabel('Wind Speed [m/s]')
plt.ylabel('Peak Temperature [C]')

plt.show()
