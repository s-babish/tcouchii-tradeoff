##!/usr/bin/env python
## Created by Julie Allen Ph.D. and Robert del Carlo
## print out one table per animal with the values from the 7 files as below
## within each table will be information from 7 files
##
##	filename = REdC10-M1-Rheo501.dat.csv
##		animal ID = REdC10-M1
##		muscle = M1
##		protocol = Rheo
##		pulse witdth (us) = (501 - 1)/10 = 50
##
##		first check that the stimulus column does not = 0
##			if It does then use set row numbers
##			if NOT then read through stimulus column
##				find the first instance where the stimulus is greater than 2
##				then from 1 less than  position take the average of the PRIOR 100 cells of FORCE column = A
##				then starting at cell where stimulus is > 2 take 1500 cells forward and find the maximum Force in that mix = B
##				Take B - A and goes into the table
## Column headers in input data files must be identical (Time, Force, Position, Stimulus)
## Create output file as such:
##		SnakeID	Muscle	Stimulus (mA)	Pulse_Width	Max_Force(g)
##		REdC10	M1		10				50			result
##		REdC10	M1		20				50			result
##		REdC10	M1		30				50			result
##		.
##		.
##		.
##		REdC10	M1		400				50000		result
##		REdC10	M1		500				50000		result
##		REdC10	M1		600				50000		result
##### EDIT THE DIRECTORY (LINE 47) SO IT JUST USES THE CURRENT DIRECTORY WHERE DATA & PYTHON SCRIPT ARE SAVED

import csv
from decimal import *
import statistics
import re
import os

my_stimulus_print_out=[10,20,30,40,50,60,70,80,90,100,200,300,400,500,600]
my_row_number_list = [33562,66893,100224,133560,166891,200228,233558,266894,300225,333561,366891,400227,433557,466893,500229]

## FIRST STEP:  find all the animals in a folder then for each animal go throught the process and create a table
## read in a list of files for an animal find all the unique animal names
experiment_list=[]
experiment_dict={}
my_stimulus=()
my_dir='C:/Users/sdbab/Downloads/Tcouchii_side_hustle/Tcouchii_side_hustle/Couchii_Rheobase'
for filename in os.listdir(my_dir):
	if filename.endswith(".dat.csv"):
		m = re.search('((^\S+?)-(\S+?))-(\D+)(\d+)', filename)
		experiment=str(m.group(1))
		muscle=str(m.group(3))
		indicated=int(m.group(5))
		pulsewidth=int((indicated - 1)/10)
		experiment_dict[pulsewidth] = filename
		experiment_list.append(experiment)
print("all the unique experiments are: ")
unique=set(experiment_list)
print(unique)

### SECOND STEP -- for each unique experiment (snake-muscle combination) get all the files, calculate the pulse width for each experiment and order the files by the pulse width
for each_experiment in unique:
	experiment_file_list=[]
	outfile_name=each_experiment + '.rheobase.csv'
	outfile=open(outfile_name, 'w')
	outfile.write("SnakeID,Muscle,Stimulus(mA),Pulse_Width(us),Max_Force(g)\n")
	print("unique experiment :" + each_experiment)
	for filename in os.listdir(my_dir):
		if filename.endswith(".dat.csv"):
			if each_experiment in filename:
				#REdC11-M1-Rheo100001.dat.csv
				print(filename)
				m = re.search('(^\S+?)-(\S+?)-(\D+)(\d+)', filename)
				my_number=int(m.group(4))
				experiment_file_list.append(my_number)
				experiment_dict[pulsewidth] = filename
	experiment_file_list.sort()
	Formatted_PW_list=[]
	for each in experiment_file_list:
		pulsewidth = int((each-1)/10)
		Formatted_PW_list.append(pulsewidth)
	print(Formatted_PW_list)
	unique_file_list=set(experiment_file_list)
	for number in experiment_file_list:
		for filename in os.listdir(my_dir):
			if filename.endswith(".dat.csv"):
				if each_experiment in filename:
					exp_combo = each_experiment.split("-")
					snake = exp_combo[0]
					muscle_num = exp_combo[1]
					if str(number) in filename:
						print("reading through the following file: ")
						print(filename)
						m = re.search('((^\S+?)-(\S+?))-(\D+)(\d+)', filename)
						experiment=m.group(1)
						muscle=m.group(3)
						protocol=m.group(4)
						indicated=int(m.group(5))
						pulsewidth=int((indicated - 1)/10)
						with open(filename) as csvfile:
							reader=csv.DictReader(csvfile)
							sum_stim = Decimal(0)
							stim=[]
							force=[]
							for row in reader:
								sum_stim += Decimal(row['Stimulus'])
								stim.append(Decimal(row['Stimulus']))
								force.append(Decimal(row['Force (g)']))
							row_num=0
							stimulus_num=0
							if sum_stim > 0:
								line_pass = 10000
								for each in stim:
									row_num = row_num + 1
									if each > 2:
										if row_num > line_pass:
											line_pass = row_num
											location = row_num - 1
											section = location - 1000
											baseline_force_list = force[section:location]
											average_baseline_force = (sum(baseline_force_list)/1000)
											max_force_list = force[(location + 1):(location + 1501)]
											maximum_force = max(max_force_list)
											force_amplitude = maximum_force - average_baseline_force #Maximum force is reported in (g), uncorrected for muscle mass
											stimnum = my_stimulus_print_out[stimulus_num]
											my_print_list = [snake,muscle_num,stimnum,pulsewidth,force_amplitude]
											for each in my_print_list:
												outfile.write(str(each))
												outfile.write(',')
											print(snake,muscle_num,stimnum,pulsewidth,force_amplitude)
											outfile.write('\n')
											line_pass = line_pass + 10000
											stimulus_num=stimulus_num+1
							elif sum_stim == 0:
								list_pos=0
								for each in stim:
									row_num = row_num + 1
									list_len = len(my_row_number_list)
									if list_pos < list_len:
										if row_num == my_row_number_list[list_pos]:
											location = row_num - 1
											section = location -1000
											baseline_force_list = force[section:location]
											average_baseline_force = (sum(baseline_force_list)/1000)
											max_force_list = force[(location+1):(location + 1501)]
											maximum_force = max(max_force_list)
											force_amplitude = maximum_force - average_baseline_force #Reported in (g), uncorrected for mass
											stimnum=my_stimulus_print_out[stimulus_num]
											print(snake,muscle_num,stimnum,pulsewidth,force_amplitude)
											my_print_list = [snake,muscle_num,stimnum,pulsewidth,force_amplitude]
											for each in my_print_list:
												outfile.write(str(each))
												outfile.write(',')
											outfile.write('\n')
											list_pos = list_pos + 1
											stimulus_num=stimulus_num+1
