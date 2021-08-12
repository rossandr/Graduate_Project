import numpy as np
import matplotlib.pyplot as plt
import xlrd

file_location = ('G:/Thermocouple Data/Ross Downward/Ross_Smolder_Test1_20201218.lvm.xls')
wb = xlrd.open_workbook(file_location)
sheet = wb.sheet_by_index(0)
sheet.cell_value(0, 0)


Thermocouple10 = []
Thermocouple11 = []
Thermocouple13 = []
Thermocouple14 = []
Thermocouple15 = []
Thermocouple16 = []
Thermocouple17 = []
Thermocouple18 = []
time = []

## Horizontal Tests ##
# Test_number = 7
# Test_number = Test_number - 1
#
# long = np.asarray([5269, 4000, 5895, 4365, 4525, 3936, 4905, 4108, 4708, 4225, 3437, 5053, 3103, 2998])
# How_long = long[Test_number]
# length = range(22, How_long)  #Change this value for each Test (keep, change)

### Downward Tests ###
Test_number = 1
Test_number = Test_number - 1

long = np.asarray([8632])
How_long = long[Test_number]
length = range(22, How_long)  #Change this value for each Test (keep, change)

def read_excell(row2, list, length):

    for step in length:
        value = sheet.cell_value(step, row2)
        list.append(value)

    return list

### Horizontal Tests small reactor###
Thermocouple10 = read_excell(21, Thermocouple10, length)  #TC 10
Thermocouple11 = read_excell(23, Thermocouple11, length)  # TC 11
Thermocouple13 = read_excell(27, Thermocouple13, length)  # TC 13
Thermocouple14 = read_excell(29, Thermocouple14, length)  # TC 14
Thermocouple15 = read_excell(31, Thermocouple15, length)  # TC 15
Thermocouple16 = read_excell(33, Thermocouple16, length)  # TC 16
Thermocouple17 = read_excell(35, Thermocouple17, length)  # TC 17
Thermocouple18 = read_excell(37, Thermocouple18, length)  # TC 17

for i in length:
    t_value = sheet.cell_value(i, 20)
    time.append(t_value)

time = np.asarray(time)

################# Plots #########################

plt.figure(1, figsize=(7, 7))
plt.plot((time)/60, Thermocouple10, label='TC 1')
plt.plot((time)/60, Thermocouple11, label='TC 2')
plt.plot((time)/60, Thermocouple13, label='TC 3')
plt.plot((time)/60, Thermocouple14, label='TC 4')
plt.plot((time)/60, Thermocouple15, label='TC 5')
plt.plot((time)/60, Thermocouple16, label='TC 6')
plt.plot((time)/60, Thermocouple17, label='TC 7')
plt.plot((time)/60, Thermocouple18, label='TC 8')
plt.xlabel('Time [min]')
plt.ylabel('Temperature [C]')
plt.grid()
plt.xlim(0)
plt.legend()

plt.figure(2, figsize=(7, 7))
plt.plot(time, Thermocouple10, label='TC 1')
plt.plot(time, Thermocouple11, label='TC 2')
plt.plot(time, Thermocouple13, label='TC 3')
plt.plot(time, Thermocouple14, label='TC 4')
plt.plot(time, Thermocouple15, label='TC 5')
plt.plot(time, Thermocouple16, label='TC 6')
plt.plot(time, Thermocouple17, label='TC 7')
plt.plot(time, Thermocouple18, label='TC 8')
plt.xlabel('Time [s]')
plt.ylabel('Temperature [C]')
plt.grid()
plt.xlim(0)
plt.legend()
# plt.show()

##################### Peak Temperatures ###########################

Temperatures = np.asarray([Thermocouple10, Thermocouple11, Thermocouple13, Thermocouple14, Thermocouple15, Thermocouple16, Thermocouple17, Thermocouple18])
peak_Temps = []

for i in range(len(Temperatures)):
    peak_Temps.append(np.max(Temperatures[i]))

print('Peak Temperatures')
print(peak_Temps)
print('Max Temperature')
print(np.max(peak_Temps))
index = np.where(peak_Temps == np.max(peak_Temps))
index = index[0] + 8
print('Max temp is %(temp)d at themocouple %(num)d' % {'temp':np.max(Temperatures), 'num': index})


################################ Spread Rates of TC 9 to 14########################################

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

set_temp = 300  #Set temperature

Thermocouple11_arr = np.asarray(Thermocouple11)
Thermocouple13_arr = np.asarray(Thermocouple13)
Thermocouple14_arr = np.asarray(Thermocouple14)
Thermocouple15_arr = np.asarray(Thermocouple15)
Thermocouple16_arr = np.asarray(Thermocouple16)
Thermocouple17_arr = np.asarray(Thermocouple17)
Thermocouple18_arr = np.asarray(Thermocouple18)


TC11_value, TC11_time_index = find_nearest(Thermocouple11_arr, set_temp)
TC14_value, TC14_time_index = find_nearest(Thermocouple14_arr, set_temp)
TC15_value, TC15_time_index = find_nearest(Thermocouple15_arr, set_temp)
TC18_value, TC18_time_index = find_nearest(Thermocouple18_arr, set_temp)
TC17_value, TC17_time_index = find_nearest(Thermocouple17_arr, set_temp)

print('\nCheck time [seconds]')
print('time index 2')
print(TC11_time_index)
print('time index 5')
print(TC15_time_index)
print('time index 7')
print(TC17_time_index)
print('time index 8')
print(TC18_time_index)

dist_11_18 = 60 # mm
dist_11_17 = 50
dist_11_15 = 30
dist_11_14 = 20
# dt = (time[TC18_time_index] - time[TC11_time_index])/60
# dt = (time[TC17_time_index] - time[TC11_time_index])/60
dt = (time[TC14_time_index] - time[TC11_time_index])/60


spread_rate = dist_11_14/dt
print('\nSpread rate mm/min')
print(spread_rate)


plt.show()




