from __future__ import print_function
import subprocess
import os
import sys
import math
import matplotlib.pyplot as plt
import numpy
import seaborn as sns; sns.set(color_codes=True)




spacer = '   '
print(os.linesep*2)
print(spacer + '____       __  ___ ___ ____ ____ ')
print(spacer + '|__/ |    |  |  |   |  |___ |__/ ')
print(spacer + '|    |___ |__|  |   |  |___ |  \ ')
print(os.linesep*2)

#Define parameters.
min_side_size = 5
max_side_size = 6
min_time = 2
max_time = 2
version_list = ['asm', 'cpp', 'asm_omp', 'cpp_omp', 'icc']
test_folder = '/home/martin/Desktop/repos/navier-stokes-simd-solver/2D/tests/to_plot/'







def read_as_array(min_time, max_time, min_side_size, max_side_size, version):



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

	data_array=(numpy.array([numpy.array(xi) for xi in data])).transpose()
	return data_array

result_asm = read_as_array(min_time, max_time, min_side_size, max_side_size, version_list[0])
#result_cpp = read_as_array(min_time, max_time, min_side_size, max_side_size, version_list[1])
#result_cpp_omp = read_as_array(min_time, max_time, min_side_size, max_side_size, version_list[3])

ax = sns.tsplot(data=result_asm, ci='sd', color='red', condition = 'ASM')
#ax2 = sns.tsplot(data=result_cpp, ci='sd', color='green', condition = 'C++')
#ax3 = sns.tsplot(data=result_cpp_omp, ci='sd', color='blue', condition = 'OpenMP')



plt.savefig('foo.png', dpi=900)


plt.show()