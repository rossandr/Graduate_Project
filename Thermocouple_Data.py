import numpy as np
import matplotlib.pyplot as plt
import xlrd

file_location = ('G:/Thermocouple Data/ME552 Ross/ME552Project_Test7_cellulose_Ross_txt.xls')
wb = xlrd.open_workbook(file_location)
sheet = wb.sheet_by_index(0)
sheet.cell_value(0, 0)


Thermocouple8 = []
Thermocouple9 = []
Thermocouple10 = []
Thermocouple11 = []
Thermocouple12 = []
Thermocouple13 = []
Thermocouple14 = []
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

long = np.asarray([19727])
How_long = long[Test_number]
length = range(22, How_long)  #Change this value for each Test (keep, change)

def read_excell(row2, list, length):

    for step in length:
        value = sheet.cell_value(step, row2)
        list.append(value)

    return list

### Horizontal Tests ###
Thermocouple8 = read_excell(17, Thermocouple8, length)
Thermocouple9 = read_excell(19, Thermocouple9, length)
Thermocouple10 = read_excell(21, Thermocouple10, length)
Thermocouple11 = read_excell(23, Thermocouple11, length)
Thermocouple12 = read_excell(25, Thermocouple12, length)
Thermocouple13 = read_excell(27, Thermocouple13, length)
Thermocouple14 = read_excell(29, Thermocouple14, length)

for i in length:
    t_value = sheet.cell_value(i, 16)
    time.append(t_value)

time = np.asarray(time)

################# Plots #########################

plt.figure(1, figsize=(7, 7))
plt.plot((time)/60, Thermocouple8, label='TC 1')
plt.plot((time)/60, Thermocouple9, label='TC 2')
plt.plot((time)/60, Thermocouple10, label='TC 3')
plt.plot((time)/60, Thermocouple11, label='TC 4')
plt.plot((time)/60, Thermocouple12, label='TC 5')
plt.plot((time)/60, Thermocouple13, label='TC 6')
plt.plot((time)/60, Thermocouple14, label='TC 7')
plt.xlabel('Time [min]')
plt.ylabel('Temperature [C]')
plt.grid()
plt.xlim(0)
plt.legend()

plt.figure(2, figsize=(7, 7))
plt.plot(time, Thermocouple8, label='TC 1')
plt.plot(time, Thermocouple9, label='TC 2')
plt.plot(time, Thermocouple10, label='TC 3')
plt.plot(time, Thermocouple11, label='TC 4')
plt.plot(time, Thermocouple12, label='TC 5')
plt.plot(time, Thermocouple13, label='TC 6')
plt.plot(time, Thermocouple14, label='TC 7')
plt.xlabel('Time [s]')
plt.ylabel('Temperature [C]')
plt.grid()
plt.xlim(0)
plt.legend()
# plt.show()

##################### Peak Temperatures ###########################

Temperatures = np.asarray([Thermocouple8, Thermocouple9, Thermocouple10, Thermocouple11, Thermocouple12, Thermocouple13, Thermocouple14])
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

Thermocouple9_arr = np.asarray(Thermocouple9)
Thermocouple10_arr = np.asarray(Thermocouple10)
Thermocouple11_arr = np.asarray(Thermocouple11)
Thermocouple12_arr = np.asarray(Thermocouple12)
Thermocouple13_arr = np.asarray(Thermocouple13)
Thermocouple14_arr = np.asarray(Thermocouple14)


TC9_value, TC9_time_index = find_nearest(Thermocouple9_arr, set_temp)
TC14_value, TC14_time_index = find_nearest(Thermocouple14_arr, set_temp)
TC13_value, TC13_time_index = find_nearest(Thermocouple13_arr, set_temp)

print('\nCheck time [seconds]')
print(TC9_time_index)
print(TC13_time_index)
print(TC14_time_index)

dist_9_14 = 85 # mm
dist_9_13 = 70
dt = (time[TC14_time_index] - time[TC9_time_index])/60
# dt = (time[TC13_time_index] - time[1501])/60

spread_rate = dist_9_14/dt
print('\nSpread rate mm/min')
print(spread_rate)

######################## Spread rate between Individual TC ################################

# Temperatures = np.asarray([Thermocouple8, Thermocouple9, Thermocouple10, Thermocouple11, Thermocouple12, Thermocouple13, Thermocouple14])
#
# dist_9_10 = 21
# dist_10_11 = 17
# dist_11_12 = 16
# dist_12_13 = 16
# dist_13_14 = 14
#
# distance_array = np.asarray([dist_9_10, dist_10_11, dist_11_12, dist_12_13, dist_13_14])
#
# def spread_rate(distance_array, set_temp, Temperatures, time):
#     TC9_value, TC9_time_index = find_nearest(Temperatures[0], set_temp)
#     TC10_value, TC10_time_index = find_nearest(Temperatures[1], set_temp)
#     TC11_value, TC11_time_index = find_nearest(Temperatures[2], set_temp)
#     TC12_value, TC12_time_index = find_nearest(Temperatures[3], set_temp)
#     TC13_value, TC13_time_index = find_nearest(Temperatures[4], set_temp)
#     TC14_value, TC14_time_index = find_nearest(Temperatures[5], set_temp)
#
#     time_index = np.asarray([TC9_time_index, TC10_time_index, TC11_time_index, TC12_time_index, TC13_time_index, TC14_time_index])
#
#     spread_rate_list = []
#     for i in range(len(distance_array)):
#
#         dt = (time[time_index[i+1]] - time[time_index[i]]) / 60
#         spread_rt = distance_array[i]/dt
#         spread_rate_list.append(spread_rt)
#
#     return spread_rate_list
#
# spread_rate_list = spread_rate(distance_array, set_temp, Temperatures, time)
# print('\nSpread rate List [mm/min]')
# print(spread_rate_list)

plt.show()




