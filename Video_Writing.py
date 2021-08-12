import cv2
import numpy as np
import glob

img_array = []

for filename in glob.glob('G:/Cellulose Burn Images/Cellulose 3ms test 9/*.jpg'):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width, height)
    img_array.append(img)

##################### Writing Video at 30 FPS ###############################################

out1 = cv2.VideoWriter('Cellulose_test9_3ms_30.mp4v', cv2.VideoWriter_fourcc(*'mp4v'), 30, size)
out2 = cv2.VideoWriter('Cellulose_test9_3ms_15.mp4v', cv2.VideoWriter_fourcc(*'mp4v'), 15, size)


for i in range(len(img_array)):
    out1.write(img_array[i])
    out2.write(img_array[i])

out1.release()
out2.release()

