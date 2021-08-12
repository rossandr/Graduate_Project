
Measure_x = 17.5
Measure_y = 16.5
Measure_z = 1 + (0/8)
Weight_sample = 1276.4
Sample_MC = 27.845


Box_weight = 633.5
Weight_sample = Weight_sample - Box_weight

### For corrected bulk desnity for moisture contnet ####

Sample_MC = Sample_MC / 100
Sample_MC = 1 - Sample_MC

Weight_sample2 = Weight_sample * Sample_MC
print(Weight_sample2)

Measure_x = Measure_x * 0.0254
Measure_y = Measure_y * 0.0254
Measure_z = Measure_z * 0.0254




Weight_sample = Weight_sample * 0.001
Weight_sample2 = Weight_sample2 * 0.001



Dimension = Measure_x * Measure_y * Measure_z
Weight = Weight_sample
Weight2 = Weight_sample2

PD = Weight/Dimension
PD2 = Weight2/Dimension

print('The Packing Density is %f [kg/m^3]' %PD)
print()
print('The Corrected MC Packing Density is %f [kg/m^3]' %PD2)
print()
print('The Duff thickness is %f [cm]' % (Measure_z*100))
print()

