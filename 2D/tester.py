from __future__ import print_function
import subprocess
import os
import sys

param_file_name = "parameters.txt"
min_side_size = 2
max_side_size = 2
min_time = 10
max_time = 30
repetitions = 100

#all_versions:
#'asm','asm_o1','asm_o2','asm_o3','asm_ofast', 
#'cpp','cpp_o1','cpp_o2','cpp_o3','cpp_ofast',
#'cpp', 'asm_omp', 'cpp_omp',
#'icc', 'icc_o1', 'icc_o2', 'icc_o3', 'icc_ofast'
# 'cpp_omp_o1', 'cpp_omp_o2', 'cpp_omp_o3', 'cpp_omp_ofast'

version_list = ['cpp_omp_ofast']
#Falta desde cpp_omp_ofast. Revisar los que dieron mal.
spacer = '   '

print(spacer + '___ ____  ___ ___ ____ ___  ')
print(spacer + ' |  |___ (__   |  |___ |__) ')
print(spacer + ' |  |___ ___)  |  |___ |  \ ')
                            
total_tests = (max_time - min_time +1)*(max_side_size - min_side_size +1)*len(version_list)
current_test = 0
for time in xrange(min_time, max_time+1):
	for side_size in xrange(min_side_size, max_side_size+1):
		for version in version_list:
			param_file = open(param_file_name,"w") 

			param_file.write("xMax " + str(side_size) + os.linesep)
			param_file.write("yMax " + str(side_size) + os.linesep)
			param_file.write("tMax " + str(time) + os.linesep)
			param_file.write("nu 0.01 " + os.linesep)
			param_file.write("rho 1.0 " + os.linesep)
			param_file.write("dx 0.05" + os.linesep)
			param_file.write("dy 0.05" + os.linesep)
			param_file.write("dt 0.005" + os.linesep)
			param_file.write("al 0.5" + os.linesep)
			param_file.write("#printPercentageSteps 20" + os.linesep)
			param_file.write("percentageStop 100.0" + os.linesep)
			param_file.write("F 10.0" + os.linesep)
			param_file.write("#printWork false" + os.linesep)
			param_file.write("C_d 2.0" + os.linesep)
			param_file.write("#printPercentage true" + os.linesep)
			param_file.write("only#printFan false" + os.linesep)
			param_file.write("plotPressureInstead false" + os.linesep)

			param_file.close()

			call = '/usr/bin/time --output ./tests/time_' + version + '_size_' + str(side_size) + '_time_' + str(time) + '.txt --append -f "%e" ./main_' + version

			current_test += 1
			print(os.linesep + spacer + "Starting test set:")
			print(2*spacer + "Test " + str(current_test) + " of " + str(total_tests))
			print(2*spacer + "side_size: " + str(side_size) + " of " + str(max_side_size))
			print(2*spacer + "t_max: " + str(time))
			print(2*spacer + "repetitions: " + str(repetitions))
			print(2*spacer + "call: " + call)

			devnull = open(os.devnull, 'w')
			for r in xrange(0,repetitions):
				print(2*spacer + 'status: ' + str(float(r)*100.0/float(repetitions)) + "%" + spacer*3, end='\r')
				sys.stdout.flush()
				try:
					subprocess.call(call, shell=True,  stdout=devnull, stderr=devnull)
				except:
					print("Warning: Either user asked to end or there was an error in subprocess call.")
					exit()

			print(2*spacer + 'status: completed' + spacer*3)
			devnull.close()



