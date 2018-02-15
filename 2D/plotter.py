from __future__ import print_function
import subprocess
import os
import sys
import math


spacer = '   '
print(os.linesep*2)
print(spacer + '____       __  ___ ___ ____ ____ ')
print(spacer + '|__/ |    |  |  |   |  |___ |__/ ')
print(spacer + '|    |___ |__|  |   |  |___ |  \ ')
print(os.linesep*2)

#Define parameters.
min_side_size = 1
max_side_size = 17
min_time = 2
max_time = 2
version_list = ['asm', 'cpp', 'asm_omp', 'cpp_omp', 'icc']
version = version_list[3]
test_folder = '/home/martin/Desktop/repos/navier-stokes-simd-solver/2D/tests/cpp_omp_times_no_output_20/'
total_tests = (max_time - min_time +1)*(max_side_size - min_side_size +1)


#Fetch data files according to parameters.
data = [[] for i in range(0, total_tests)]
current_test = 0
for time in xrange(min_time, max_time+1):
	for side_size in xrange(min_side_size, max_side_size+1):
		filename = 	'time_' + version + '_size_' + str(side_size) + '_time_' + str(time) + '.txt'
		file = test_folder + filename
		file_data = []
		with open(file) as fp:
			for line in fp:
				file_data.append(float(line))

		for m in file_data:
			data[current_test].append(m)
		current_test += 1

#Calculate means
means = []
for test in data:
	means.append(float(sum(test))/float(len(test)))

#Calculate variances
variances = []
for data_idx in xrange(0,len(data)):
	sqdifs = []
	var = 0
	for elem in data[data_idx]:
		var += (means[data_idx] - elem)**2
	variances.append(var)

std_devs = [math.sqrt(var) for var in variances]
		
print(spacer + 'Means: ' + os.linesep + spacer + str(means))
print(spacer + 'Variances: ' + os.linesep + spacer + str(variances))
print(spacer + 'Stardard deviations: ' + os.linesep + spacer + str(std_devs))


import matplotlib.pyplot as plt
import numpy
import seaborn as sns; sns.set(color_codes=True)





data_array=numpy.array([numpy.array(xi) for xi in data])
print(data_array)


ax = sns.tsplot(data=data_array.transpose(), ci=[100])


plt.show()