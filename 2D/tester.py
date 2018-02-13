from __future__ import print_function
import subprocess
import os
import sys

param_file_name = "parameters.txt"
min_side_size = 1
max_side_size = 100
t_max = 2
repetitions = 105



print('  ___ ____ ____ ___ ____ ____ ')
print('   |  |___ [__   |  |___ |__/ ')
print('   |  |___ ___]  |  |___ |  \ ')
                            


for side_size in xrange(min_side_size, max_side_size+1):
	
	param_file = open(param_file_name,"w") 

	param_file.write("xMax " + str(side_size) + ".0" + os.linesep)
	param_file.write("yMax " + str(side_size) + ".0" + os.linesep)
	param_file.write("tMax 2.0" + os.linesep)
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

	call = '/usr/bin/time --output ./time_cppomp_size_' + str(side_size) + '_time_' + str(t_max) + '.txt --append -f "%e" ./mainCPP_omp'
	spacer = '   '
	print(os.linesep + spacer + "Starting test set:")
	print(2*spacer + "side_size: " + str(side_size) + " of " + str(max_side_size))
	print(2*spacer + "t_max: " + str(t_max))
	print(2*spacer + "repetitions: " + str(repetitions))
	print(2*spacer + "call: " + call)


	devnull = open(os.devnull, 'w')
	for x in xrange(1,repetitions):
		print(2*spacer + 'percentage: ' + str(float(x)*100.0/float(repetitions)) + "%", end='\r')
		sys.stdout.flush()
		subprocess.call(call, shell=True,  stdout=devnull, stderr=devnull)
	print(2*spacer + 'percentage: 100%')
	devnull.close()



