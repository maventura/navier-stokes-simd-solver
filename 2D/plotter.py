from __future__ import print_function
import subprocess
import os
import sys
import math
import matplotlib.pyplot as plt
import numpy
import seaborn as sns; sns.set(color_codes=True)
from collections import OrderedDict



spacer = '   '
print(os.linesep*2)
print(spacer + '____       __  ___ ___ ____ ____ ')
print(spacer + '|__/ |    |  |  |   |  |___ |__/ ')
print(spacer + '|    |___ |__|  |   |  |___ |  \ ')
print(os.linesep*2)

#Define parameters.
min_side_size = 2
max_side_size = 2
min_time = 1
max_time = 30

#Executable names:
version_list_asm = ['asm','asm_o3','asm_ofast']
version_list_cpp = ['cpp','cpp_o1','cpp_o2','cpp_o3','cpp_ofast']
version_list_icc = ['icc', 'icc_o1', 'icc_o2', 'icc_o3', 'icc_ofast']
version_list_omp = ['cpp', 'cpp_omp', 'cpp_omp_o1', 'cpp_omp_o2', 'cpp_omp_o3', 'cpp_omp_ofast']
version_list_ofast = ['asm_ofast','cpp_ofast','icc_ofast','cpp_omp_ofast']
version_list_all = version_list_asm+version_list_cpp+version_list_icc+version_list_omp

#Corresponding plot labels:
version_names_list_asm = ['Assembler (NASM) and GCC C++','Assembler (NASM) and GCC C++ -O3','Assembler (NASM) and GCC C++ -Ofast']
version_names_list_cpp = ['GCC C++','GCC C++ -O1','GCC C++ -O2','GCC C++ -O3','GCC C++ -Ofast']
version_names_list_icc = ['Intel C++ Compiler', 'Intel C++ Compiler -O1', 'Intel C++ Compiler -O2', 'Intel C++ Compiler -O3', 'Intel C++ Compiler -Ofast']
version_names_list_omp = ['GCC C++', 'GCC C++ OpenMP', 'GCC C++ OpenMP -O1', 'GCC C++ OpenMP -O2', 'GCC C++ OpenMP -O3', 'GCC C++ OpenMP -Ofast']
version_names_list_ofast = ['Assembler (NASM) and GCC C++ -Ofast', 'GCC C++ -Ofast', 'Intel C++ Compiler -Ofast', 'GCC C++ OpenMP -Ofast']
version_names_list_all = version_names_list_asm+version_names_list_cpp+version_names_list_icc+version_names_list_omp

#delete repeted elements from the 'all' lists preserving order
version_list_all = list(OrderedDict.fromkeys(version_list_all))
version_names_list_all = list(OrderedDict.fromkeys(version_names_list_all))




#version_list = version_list_asm
#version_names_list = version_names_list_asm

#version_list = version_list_cpp
#version_names_list = version_names_list_cpp

#version_list = version_list_icc
#version_names_list = version_names_list_icc

#version_list = version_list_omp
#version_names_list = version_names_list_omp

version_list = version_list_ofast
version_names_list = version_names_list_ofast

#version_list = version_list_all
#version_names_list = version_names_list_all

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
					try:
						file_data.append(float(line))
					except Exception as e:
						print("Error parsing file " + str(file))
						exit()

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

#plt.switch_backend('GTK')

#[u'pgf', u'ps', u'Qt4Agg', u'GTK', u'GTKAgg',
#u'nbAgg', u'agg', u'cairo', u'MacOSX', u'GTKCairo',
#u'Qt5Agg', u'template', u'WXAgg', u'TkAgg', u'GTK3Cairo',
#u'GTK3Agg', u'svg', u'WebAgg', u'pdf', u'gdk', u'WX']



sns.set_context(font_scale=0.5, rc={"lines.linewidth": 1})

for i in xrange(0,len(version_list)):
	result = read_as_array(min_time, max_time, min_side_size, max_side_size, version_list[i])
	ax = sns.tsplot(data=result, ci='sd',  color=sns.color_palette()[i%len(sns.color_palette())], condition = version_names_list[i])

plt.savefig('foo.png', dpi=900, pad_inches=0)

manager = plt.get_current_fig_manager()


plt.savefig('foo.png',  bbox_inches='tight')
#Warning: figure gets saved with a smaller frame due to a bug of matplotlib. 
#save the figure manualy from the plot window to get a better label positioning.

plt.show()

