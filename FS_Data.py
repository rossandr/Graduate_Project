from PIL import Image
import numpy as np
from skimage import io
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
import scipy as sp
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline


def calibration_image(snap_filename):
    """Used to find calibration factor of calibration image,
    The calibration factor is needed to determine the distance each pixel occupies in mm/pixel
    input: snap_filename for the burn
    returns: the coordiantes of the points and the calibration factor"""

    file_path = "F:/Smoldering Duff 2021/Tiff/IR_data_tiffs_W/"
    file = (file_path + '/' + snap_filename)
    print(type(snap_filename))

    plt.figure()
    plt.imshow(mpimg.imread(file))
    print('Click 2 times. at corner A and then corner B')
    # distance from AB and midpoint AB will use reference 0, and continue through 3 for DA
    x = plt.ginput(2)
    print(x)
    plt.close('all')

    calibration_factor = abs(200 / (np.sqrt(
        (x[0][0] - x[1][0]) ** 2 + (x[0][1] - x[1][1]) ** 2)))  # 200 is mm distance between the two plate points
    print('Calibration Factor is mm/pixel is {}'.format(calibration_factor))

    return x, calibration_factor


def cropping_for_spread_rates(tif_filename, start_frame, end_frame, x_crop_left, x_crop_right, y_crop_bottom,
                              y_crop_top):
    """Used to crop the smoldering tiff file to the desired size and length
    inputs: the tiff smoldering file, frame to start viewing, frame to end viewing, the x crop and y crop
    values found from FLIR software : the x_crop(0 - 640) y_crop(0 - 512)
    Returns im1 the croped np_array"""

    file_path = "F:/Smoldering Duff 2021/Tiff/IR_data_tiffs_W/"
    im1 = io.imread(file_path + '/' + tif_filename)

    total_frames = len(im1)
    # print('The total number of frames is: {}'.format(total_frames))

    total_x_pixel = 640
    total_y_pixel = 512
    frames_being_viewed = end_frame - start_frame
    # print('The total number of frames used is: {}'.format(frames_being_viewed))

    x_crop = (x_crop_left, x_crop_right)
    # print('The cropped x axis is: {}'.format(x_crop))

    y_crop = (y_crop_top, y_crop_bottom)
    # print('The cropped y axis is: {}'.format(y_crop))

    im1_new = im1[start_frame:end_frame, y_crop_top:y_crop_bottom, x_crop_left:x_crop_right]
    # print('The shape of the original tiff is: {}'.format(im1.shape))
    # print('The shape of the cropped tif is: {}'.format(im1_new.shape))

    return im1_new


def spread_rate(im1_new, calibration_factor=1.0196):
    """used to calculate the spread rate, assumes that the burn is circular
    inputs: the cropped np_array, the calibration factor
    returns: the spread rate (velocity), mean spread rate, mean spread rate list, and the Radius"""
    threshold = 1200
    for i in range(len(im1_new)):
        if i > 0:
            im1_new[i][im1_new[
                           i - 1] == 1] = threshold + 1  # if the previous frames photon count is above the threshold automatically sets to threshold
            im1_new[i][im1_new[i] < threshold] = 0  # if the photon count is below threshold set value to zero
            im1_new[i][im1_new[i] >= threshold] = 1  # if the photon count is above threshold set value to 1
        if i == 0:
            im1_new[i][im1_new[i] < threshold] = 0  # if the photon count is below threshold set value to zero
            im1_new[i][im1_new[i] >= threshold] = 1  # if the photon count is above threshold set value to 1

    smoldering_pixel_lst = []
    for frame in range(len(im1_new)):  # determines the number of pixels in the frame that are considered smoldering
        number_smold = np.count_nonzero(im1_new[frame] == 1.0)  # number of pixels that are == 1
        smoldering_pixel_lst.append(number_smold)
    smoldering_pixel_array = np.asarray(smoldering_pixel_lst)

    Area = smoldering_pixel_array * (
                calibration_factor ** 2)  # calculates the area that is smoldering, assuming circle, and area of square pixel
    Radius = np.sqrt(Area / np.pi)  # calculates the radius of the circle to find spread rate

    velocity = []
    dt = 10 / 60  # min
    for vel in range(len(Area) - 1):  # determines the spread velocity between consecutive frames
        velocity_ind = (Radius[vel + 1] - Radius[vel]) / dt
        velocity.append(velocity_ind)

    mean_velocity = np.mean(velocity)  # mean spread rate of all frames calculated
    mean_velocity_list = [mean_velocity] * len(velocity)
    print('the mean spread rate is {} mm/min'.format(mean_velocity))

    return velocity, mean_velocity, mean_velocity_list, Radius


def plotting_spread_rate(im1_new, calibration_factor=1.0196):
    """plots the spread rate to help reduce error, if there is a curve in the begining or end of the radius plot cut thise times out
    input: trimmed smoldering array
    returns: plot of spread rate vs time"""
    velocity, mean_velocity, mean_velocity_list, Radius = spread_rate(im1_new,
                                                                      calibration_factor)  # finds velocities form previous function

    dt = 10 / 60

    time = np.arange(0, len(velocity) * dt, dt)  # plots velocity vs time
    plt.figure(1)
    plt.plot(time, velocity, color='blue')
    plt.plot(time, mean_velocity_list, linestyle='dashed', color='red')
    plt.xlabel('Time min')
    plt.ylabel('Spread Rate mm/min')
    plt.title('Spread Rate')

    time2 = np.arange(0, len(Radius) * dt, dt)  # plots radius vs time
    plt.figure(2)
    plt.plot(time2, Radius)
    plt.xlabel('Time min')
    plt.ylabel('Radius')
    plt.title('Radius vs Time')
    plt.show()


# def intensity_to_temperature(filename):
#     """converts photon counts to temperatures
#     input: tiff smoldering file
#     returns: np array of temperatures"""
#
#     try: #checks to see if this has already been done and file has been saveed to reduce time
#         filename1 = filename.split('.')
#         filename2 = filename1[0]
#         filename3 = ('Temps_' + filename2 + '.npy')
#         temperatures = np.load(filename3)
#         print('opened existing file')
#
#     except: #if file is not found then it calculates temperatures from intensity and saves to np array to npy file
#
#         file_path = "G:/Ross_Duff_Smolder_TIFF/FS_Duff_Smolder/"
#         intensity = io.imread(file_path + '/' + filename)
#
#         # below is a black body calibration curve used to convert photon counts to temperatures
#         # this was found experimentally using the exact setup to ensure temps are accurate
#         emissivity = 0.85
#         bb_temp = np.arange(50, 650, 50)
#         bb_intensity = [828, 926, 1185, 1648, 2404, 3510, 5166, 7270, 10029, 13695, 16273, 16294]
#         bb_temp_ary1 = np.asarray(bb_temp)
#         bb_temp_ary1_K = bb_temp_ary1 + 273.15
#         bb_intensity_ary = np.asarray(bb_intensity)
#         bb_temp_ary = (bb_temp_ary1_K * 0.98) / emissivity
#         bb_temp_ary = bb_temp_ary - 273.15
#
#         linear_function = interpolate.interp1d(bb_intensity, bb_temp_ary) #interpolates BB profile to find temperatures
#         temperatures = linear_function(intensity)
#         filename1 = filename.split('.')
#         filename2 = filename1[0]
#         filename3 = ('Temps_' + filename2)
#         filename4 = ('Temps_' + filename2 + '.txt')
#         np.save(filename3, temperatures)
#         print('created new file named:' + filename3)
#         np.savetxt(filename4, temperatures[0]) #saves temperature np_array to file to reduce program cost
#
#     return temperatures

def cropping_for_intensity(tif_filename, start_frame_temps, end_frame_temps, x_crop_left_temps, x_crop_right_temps,
                       y_crop_bottom_temps, y_crop_top_temps):
    """used to crop the temperatures to the desired times and x and y frame locations
    inputs: the tiff smoldering file, frame to start viewing, frame to end viewing, the x crop and y crop
    values found from FLIR software : the x_crop(0 - 640) y_crop(0 - 512)
    Returns im1 the croped np_array """

    file_path = "F:/Smoldering Duff 2021/Tiff/IR_data_tiffs_W/"
    intensity_ary = io.imread(file_path + '/' + tif_filename)

    # intensity_ary = intensity_to_temperature(tif_filename)  # gets temps from photon count data

    intensity_new = intensity_ary[start_frame_temps:end_frame_temps, y_crop_top_temps:y_crop_bottom_temps, x_crop_left_temps:x_crop_right_temps]  # new cropped array

    return intensity_new

def intensity_to_temperature(tif_filename, start_frame_temps, end_frame_temps, x_crop_left_temps, x_crop_right_temps,
                       y_crop_bottom_temps, y_crop_top_temps):
    """converts photon counts to temperatures
    input: tiff smoldering file
    returns: np array of temperatures"""

    intensity_cropped = cropping_for_intensity(tif_filename, start_frame_temps, end_frame_temps, x_crop_left_temps, x_crop_right_temps,
                       y_crop_bottom_temps, y_crop_top_temps)

    # print(np.max(intensity_cropped))
    # print(np.min(intensity_cropped))

    threshold = 16290

    for i in range(len(intensity_cropped)):
        intensity_cropped[i][intensity_cropped[i] > threshold] = 16294  # if temp is above smolder threshold set to 0


    # below is a black body calibration curve used to convert photon counts to temperatures
    # this was found experimentally using the exact setup to ensure temps are accurate
    # the following black body test was performed with an emissivity around 0.98
    emissivity = 0.85
    # bb_temp = np.arange(50, 600, 50) #[50,100,150,200,250,300,350,400,450,500,550,600]
    bb_temp = [50*(f+1) for f in range(12)]
    print(bb_temp)

    bb_intensity = [828, 926, 1185, 1648, 2404, 3510, 5166, 7270, 10029, 13695, 16273, 16294]
    bb_temp_ary1 = np.asarray(bb_temp) # Converts BB_temps into an array
    bb_temp_ary1_K = bb_temp_ary1 + 273.15 # converts from Celsius to Kelvin
    bb_intensity_ary = np.asarray(bb_intensity)
    bb_temp_ary = (bb_temp_ary1_K * 0.98) / emissivity
    bb_temp_ary = bb_temp_ary - 273.15 # converts back from K to C

    #print('linear interpolation failed')

    linear_function = interpolate.interp1d(bb_intensity, bb_temp_ary)  # interpolates BB profile to find temperatures
    # linear_function = InterpolatedUnivariateSpline(bb_intensity, bb_temp_ary, k=3)  # interpolates BB profile to find temperatures
    #print('temperatures failed')
    temperatures = linear_function(intensity_cropped)



    return temperatures


# def cropping_for_temps(tif_filename, start_frame_temps, end_frame_temps, x_crop_left_temps, x_crop_right_temps,
#                        y_crop_bottom_temps, y_crop_top_temps):
#     """used to crop the temperatures to the desired times and x and y frame locations
#     inputs: the tiff smoldering file, frame to start viewing, frame to end viewing, the x crop and y crop
#     values found from FLIR software : the x_crop(0 - 640) y_crop(0 - 512)
#     Returns im1 the croped np_array """
#
#     temps_ary = intensity_to_temperature(tif_filename)  # gets temps from photon count data
#
#     frames_being_viewed = end_frame_temps - start_frame_temps
#     # print('The total number of frames used is: {}'.format(frames_being_viewed))
#
#     x_crop = (x_crop_left_temps, x_crop_right_temps)
#     # print('The cropped x axis is: {}'.format(x_crop))
#
#     y_crop = y_crop_bottom_temps - y_crop_top_temps
#     # print('The cropped y axis is: {}'.format(y_crop))
#
#     temps_new = temps_ary[start_frame_temps:end_frame_temps, y_crop_top_temps:y_crop_bottom_temps, x_crop_left_temps:x_crop_right_temps]  # new cropped array
#     # print('The shape of the original tiff is: {}'.format(temps_ary.shape))
#     # print('The shape of the cropped tif is: {}'.format(temps_new.shape))
#
#     return temps_new


def averaged_temperatures_max(temps_new):
    """used to calculate the average temperatures over the desired frames
    input: smolder file name
    returns: average temperatures over those frames"""
    threshold = 365

    temps_new_count = np.copy(temps_new)
    for i in range(len(temps_new_count)):
        temps_new_count[i][temps_new_count[i] < threshold] = 0.0  # if temp is below smolder threshold set to 0
        temps_new_count[i][temps_new_count[i] >= threshold] = 1.0  # if temp is above smolder threshold set to 1

    pixels_above_threshold = []
    for frames in range(len(temps_new)):  # determines how many pixels were considered smoldering
        count_above_threshold = np.count_nonzero(temps_new_count[frames] == 1.0)
        pixels_above_threshold.append(count_above_threshold)

    for i in range(len(temps_new)):
        temps_new[i][temps_new[i] < threshold] = 0

    summation = []
    for frames in range(len(temps_new)):  # sums the number of pixels considered smoldering
        summation_per_frame = np.sum(temps_new[frames])
        summation.append(summation_per_frame)

    average_temp = []
    for frames in range(len(summation)):  # finds the average temperature from desired frames and spatial distancing
        average_temp_per_frame = summation[frames] / pixels_above_threshold[frames]
        average_temp.append(average_temp_per_frame)

    return average_temp

def averaged_temperatures_bound(temps_new):
    """used to calculate the average temperatures over the desired frames
    input: smolder file name
    returns: average temperatures over those frames"""
    threshold = 220
    top_threshold = 365

    temps_new_count = np.copy(temps_new)
    for i in range(len(temps_new_count)):
        temps_new_count[i][temps_new_count[i] < threshold] = 0   # if temp is below smolder threshold set to 0
        temps_new_count[i][temps_new_count[i] > top_threshold] = 0   # if temp threshold is above 365 set to 0
        temps_new_count[i][temps_new_count[i] >= threshold] = 1  # if temp is above smolder threshold set to 1

    pixels_above_threshold = []
    for frames in range(len(temps_new)):  # determines how many pixels were considered smoldering
        count_above_threshold = np.count_nonzero(temps_new_count[frames] == 1.0)
        pixels_above_threshold.append(count_above_threshold)

    for i in range(len(temps_new)):
        temps_new[i][temps_new[i] < threshold] = 0
        temps_new[i][temps_new[i] > top_threshold] = 0

    summation = []
    for frames in range(len(temps_new)):  # sums the number of pixels considered smoldering
        summation_per_frame = np.sum(temps_new[frames])
        summation.append(summation_per_frame)

    average_temp = []
    for frames in range(len(summation)):  # finds the average temperature from desired frames and spatial distancing
        average_temp_per_frame = summation[frames] / pixels_above_threshold[frames]
        average_temp.append(average_temp_per_frame)

    return average_temp


def plotting_average_temp_bound(temps1, frame_start, frame_stop):
    """used to plot aveerage temperature over time
    inputs: np temp array
    returns: average_temps vs time plot"""
    average_temp = averaged_temperatures_bound(temps1)  # find average temp at each frame
    x_frame = list(range(frame_start, frame_stop))

    dt = 10 / 60
    time = np.arange(0, len(average_temp) * dt, dt)

    # plt.figure(3)
    # plt.plot(time, average_temp, color='red')
    # plt.xlabel('Time min')
    # plt.ylabel('Average Temperature [degrees C]')
    # plt.title('Average Temperatures Bound')


    plt.figure(4)
    plt.plot(x_frame, average_temp, color='red')
    plt.xlabel('Frame')
    plt.ylabel('Average Temperature [degrees C]')
    plt.title('Average Temperatures Bound')
    plt.show()

def plotting_average_temp_max(temps1, frame_start, frame_stop):
    """used to plot aveerage temperature over time
    inputs: np temp array
    returns: average_temps vs time plot"""
    average_temp = averaged_temperatures_max(temps1)  # find average temp at each frame
    x_frame = list(range(frame_start, frame_stop))

    dt = 10 / 60
    time = np.arange(0, len(average_temp) * dt, dt)

    # plt.figure(3)
    # plt.plot(time, average_temp, color='red')
    # plt.xlabel('Time min')
    # plt.ylabel('Average Temperature [degrees C]')
    # plt.title('Average Temperatures Max')

    plt.figure(4)
    plt.plot(x_frame, average_temp, color='red')
    plt.xlabel('Frame')
    plt.ylabel('Average Temperature [degrees C]')
    plt.title('Average Temperatures max')
    plt.show()


# def average_temp_at_frame(temps1, frame, start_frame_main):
#     """used to find average smoldering temperature at a single frame
#     inputs: temperature array, frame to consider, start croped frame
#     returns: the average temperature of smoldering at the desired frame"""
#
#     average_temp = averaged_temperatures_bound(temps1)
#     print(average_temp)
#     print(len(average_temp))
#     frame_calculated = frame - start_frame_main  # actual frame
#     print('frame calculated is {}'.format(frame_calculated))
#     average_temp_desired = average_temp[frame_calculated]
#
#     print('average temperature at frame {} is : {}'.format(frame, average_temp_desired))
#
#     ## Area burned
#     calibration_factor = 1.038669
#     area_smolder_count = np.count_nonzero(temps1[frame_calculated])
#     area = area_smolder_count * calibration_factor ** 2
#
#     print('Area at frame {} in: {} mm^2'.format(frame, area))
#
#     return average_temp_desired, area


def averaged_Watt_max(temps_new):
    """used to calculate the average temperatures over the desired frames
    input: smolder file name
    returns: average temperatures over those frames"""
    threshold = 365

    temps_new_count = np.copy(temps_new)
    for i in range(len(temps_new_count)):
        temps_new_count[i][temps_new_count[i] < threshold] = 0.0  # if temp is below smolder threshold set to 0
        temps_new_count[i][temps_new_count[i] >= threshold] = 1.0  # if temp is above smolder threshold set to 1

    pixels_above_threshold = []
    for frames in range(len(temps_new)):  # determines how many pixels were considered smoldering
        count_above_threshold = np.count_nonzero(temps_new_count[frames] == 1.0)
        pixels_above_threshold.append(count_above_threshold)

    for i in range(len(temps_new)):
        temps_new[i][temps_new[i] < threshold] = 0

    W_summation = []
    W_b_summation = []
    for frames in range(len(temps_new)):  # sums the number of pixels considered smoldering
        calibration_factor = 1.038669 * 0.001 #[m^2]
        emmissivity = 0.85
        boltzman = 5.670374419E-8
        W_b = emmissivity * boltzman * temps_new ** 4  # [Watt/m^2]
        Watt = W_b * calibration_factor**2

        W_b_summation_per_frame = np.sum(W_b[frames])
        W_summation_per_frame = np.sum(Watt[frames])
        W_summation.append(W_b_summation_per_frame)
        W_b_summation.append(W_summation_per_frame)

    average_W = []
    for frames in range(len(W_summation)):  # finds the average temperature from desired frames and spatial distancing
        average_W_per_frame = W_summation[frames] / pixels_above_threshold[frames]
        average_W.append(average_W_per_frame)

    average_W_b = []
    for frames in range(len(W_b_summation)):  # finds the average temperature from desired frames and spatial distancing
        average_W_b_per_frame = W_b_summation[frames] / pixels_above_threshold[frames]
        average_W_b.append(average_W_b_per_frame)

    return average_W, average_W_b, W_b_summation, W_summation

def averaged_Watt_bound(temps_new):
    """used to calculate the average temperatures over the desired frames
    input: smolder file name
    returns: average temperatures over those frames"""
    threshold = 220
    top_threshold = 365

    temps_new_count = np.copy(temps_new)
    for i in range(len(temps_new_count)):
        temps_new_count[i][temps_new_count[i] < threshold] = 0.0  # if temp is below smolder threshold set to 0
        temps_new_count[i][temps_new_count[i] > top_threshold] = 1.0 # if temp above top threshold consider 0
        temps_new_count[i][temps_new_count[i] >= threshold] = 1.0  # if temp is above smolder threshold set to 1

    pixels_above_threshold = []
    for frames in range(len(temps_new)):  # determines how many pixels were considered smoldering
        count_above_threshold = np.count_nonzero(temps_new_count[frames] == 1.0)
        pixels_above_threshold.append(count_above_threshold)

    for i in range(len(temps_new)):
        temps_new[i][temps_new[i] < threshold] = 0
        temps_new[i][temps_new[i] > top_threshold] = 0


    W_summation = []
    W_b_summation = []
    for frames in range(len(temps_new)):  # sums the number of pixels considered smoldering
        calibration_factor = 1.038669 * 0.001 #[m^2]
        emmissivity = 0.85
        boltzman = 5.670374419E-8
        W_b = emmissivity * boltzman * temps_new ** 4  # [Watt/m^2]
        Watt = W_b * calibration_factor**2

        W_b_summation_per_frame = np.sum(W_b[frames])
        W_summation_per_frame = np.sum(Watt[frames])
        W_summation.append(W_b_summation_per_frame)
        W_b_summation.append(W_summation_per_frame)

    average_W = []
    for frames in range(len(W_summation)):  # finds the average temperature from desired frames and spatial distancing
        average_W_per_frame = W_summation[frames] / pixels_above_threshold[frames]
        average_W.append(average_W_per_frame)

    average_W_b = []
    for frames in range(len(W_b_summation)):  # finds the average temperature from desired frames and spatial distancing
        average_W_b_per_frame = W_b_summation[frames] / pixels_above_threshold[frames]
        average_W_b.append(average_W_b_per_frame)

    return average_W, average_W_b, W_b_summation, W_summation


def radiant_emittance_max(temps1, frame, start_frame_main):
    print()
    print('################# Max Watt Data ################')
    print()

    W_average, W_b_average, W_b_summation, W_summation = averaged_Watt_max(temps1)
    frame_calculated = frame - start_frame_main  # actual frame
    average_W_desired = W_average[frame_calculated]
    print('frame calculated is {}'.format(frame_calculated))
    average_W_b_desired = W_b_average[frame_calculated]

    print('Total W at frame {} is : {} [W]'.format(frame, W_summation[frame_calculated]))
    print('average W at frame {} is : {} [W]'.format(frame, average_W_desired))
    print('Total W_b at frame {} is : {} [watt/m^2]'.format(frame, W_b_summation[frame_calculated]))
    print('average W_b at frame {} is : {} [watt/m^2]'.format(frame, average_W_b_desired))

    return

def radiant_emittance_bound(temps1, frame, start_frame_main):
    print()
    print('################# Bounded Watt Data ################')
    print()

    W_average, W_b_average, W_b_summation, W_summation = averaged_Watt_bound(temps1)
    frame_calculated = frame - start_frame_main  # actual frame
    average_W_desired = W_average[frame_calculated]
    print('frame calculated is {}'.format(frame_calculated))
    average_W_b_desired = W_b_average[frame_calculated]

    print('Total W at frame {} is : {} [W]'.format(frame, W_summation[frame_calculated]))
    print('average W at frame {} is : {} [W]'.format(frame, average_W_desired))
    print('Total W_b at frame {} is : {} [watt/m^2]'.format(frame, W_b_summation[frame_calculated]))
    print('average W_b at frame {} is : {} [watt/m^2]'.format(frame, average_W_b_desired))

    return


def average_temp_at_frame_bound(temps1, frame, start_frame_main):
    """used to find average smoldering temperature at a single frame
    inputs: temperature array, frame to consider, start croped frame
    returns: the average temperature of smoldering at the desired frame"""

    print()
    print('################# Bounded Temperature Data ################')
    print()
    average_temp = averaged_temperatures_bound(temps1)
    #print(len(average_temp))
    frame_calculated = frame - start_frame_main  # actual frame
    # print('frame calculated is {}'.format(frame_calculated))
    average_temp_desired = average_temp[frame_calculated]

    print('average temperature at frame {} is : {}'.format(frame, average_temp_desired))

    print('Max temperature is {}'.format(np.max(temps[frame_calculated])))
    print('Number of pixels considered smoldering = {}'.format(np.count_nonzero(temps1[frame_calculated])))

    initial_frame = frame - 10
    frames = list(range(initial_frame, frame1))
    print(frames)
    average_smold_multiple = []
    for i in range(len(frames)):
        frame_calculated1 = (frame - start_frame_main) - 10 + i
        average_temp_desired1 = average_temp[frame_calculated1]
        average_smold_multiple.append(average_temp_desired1)

    sum_of_temps = (sum(average_smold_multiple))/10
    print('10 frame temperature average is {}'.format(sum_of_temps))

    number_pix_smold = np.count_nonzero(temps1[frame_calculated])


    ## Area burned
    calibration_factor = 1.038669
    area_smolder_count = np.count_nonzero(temps1[frame_calculated])
    area = area_smolder_count * calibration_factor ** 2

    print('Area at frame {} in: {} mm^2'.format(frame, area))
    print('Area at frame {} in: {} cm^2'.format(frame, area * 0.01))

    return average_temp_desired, area, number_pix_smold

def average_temp_at_frame_max(temps1, frame, start_frame_main):
    """used to find average smoldering temperature at a single frame
    inputs: temperature array, frame to consider, start croped frame
    returns: the average temperature of smoldering at the desired frame"""

    print()
    print('################# Max Temperature Data ################')
    print()

    average_temp = averaged_temperatures_max(temps1)
    # print(len(average_temp))
    frame_calculated = frame - start_frame_main  # actual frame
    # print('frame calculated is {}'.format(frame_calculated))
    average_temp_desired = average_temp[frame_calculated]

    print('average temperature at frame {} is : {}'.format(frame, average_temp_desired))

    print('Max temperature is {}'.format(np.max(temps[frame_calculated])))
    print('Number of pixels considered smoldering = {}'.format(np.count_nonzero(temps1[frame_calculated])))
    number_pix_smold = np.count_nonzero(temps1[frame_calculated])

    ## Area burned
    calibration_factor = 1.038669
    area_smolder_count = np.count_nonzero(temps1[frame_calculated])
    area = area_smolder_count * calibration_factor ** 2

    print('Area at frame {} in: {} mm^2'.format(frame, area))
    print('Area at frame {} in: {} cm^2'.format(frame, area * 0.01))

    return average_temp_desired, area, number_pix_smold


if __name__ == '__main__':
    tiff = "20210510_T1_smold"
    snap = "20210510_T1_smold"

    filename_tif = (tiff + '.tif')
    filename_snap = (snap + '.tif')

    # desired start and stop frames
    #start_frame = 250
    #end_frame = 650

    # desired left and right crops (from left to right)
    #x_crop_left = 0
    #x_crop_right = 555

    # desired y crops (from top to bottom)
    #y_crop_top = 0
    #y_crop_bottom = 127

    # im1_new1 = cropping_for_spread_rates(filename_tif, start_frame, end_frame, x_crop_left, x_crop_right, y_crop_bottom,
    #                                      y_crop_top)
    #
    # # Use with calibration image first to get calibration_factor
    #x, calibration_factor = calibration_image(filename_snap)
    # #
    # # #use for calibration factor unknown
    # plotting_spread_rate(im1_new1, calibration_factor)

    # Use if calibration factor known
    # plotting_spread_rate(im1_new1)

    #####################################################################################################################
    # the following is for temperature calcs

    # desired start and stop frames
    start_frame_temps = 98
    end_frame_temps = 151

    # desired left and right crops (from left to right)
    x_crop_left_temps = 0  #0
    x_crop_right_temps = 640   #639

    # desired y crops (from top to bottom)
    y_crop_top_temps = 0  #0
    y_crop_bottom_temps = 512   #511

    temps = intensity_to_temperature(filename_tif, start_frame_temps, end_frame_temps, x_crop_left_temps,x_crop_right_temps, y_crop_bottom_temps, y_crop_top_temps)

    # # for temperature calculations


    ##### Plots ####
    # plotting_average_temp_bound(temps, start_frame_temps, end_frame_temps)
    # plotting_average_temp_max(temps, start_frame_temps, end_frame_temps)

    ##### Bounded Data ######
    frame1 = 150
    average_temp_at_frame_bound(temps, frame1, start_frame_temps)
    radiant_emittance_bound(temps, frame1, start_frame_temps)



    ###### Max Data #######

    average_temp_at_frame_max(temps, frame1, start_frame_temps)
    radiant_emittance_bound(temps, frame1, start_frame_temps)
