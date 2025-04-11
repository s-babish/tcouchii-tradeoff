							##!/usr/bin/env python
### Created by Julie Allen Ph.D. and modified by Robert del Carlo
### Print out one table per animal with with the maximum force values from 15 locations in the file
### These values will be plotted against known stimulus magnitudes (10(d10)100, 100(d100)600)
### This will be done for six files recorded from each individual; each at a different temperature
### The temperatures are: 10, 15, 20, 25, 30, and 35 (ºC) - they were performed in random order
### This procedure may have been performed at multiple stimulus pulse widths for optimization
##
##	filename = REdC10-M1-Rheo10000-35C.dat.csv
##		animal ID = REdC10-M1   ## snake-muscle combinations constitute one experiment
##		muscle = M1
##		protocol = Rheo10000
##		pulse witdth = 10000
##		Temperature = 35(ºC)
##
##		first check that the stimulus column does not = 0 , experiments of pulse width 50 or 100 microseconds will be zero
##			if It does then use set row numbers
##			if NOT then read through stimulus column for greater accuracy
##				find the first instance where the stimulus is greater than 2
##				then from 1 less than  position take the average of the PRIOR 100 cells of FORCE column = A
##				then starting at cell where stimulus is > 2 take 1500 cells forward and find the maximum Force in that mix = B
##				Report B - A (maximum force) into the table
##
##  Report one table per animal with all maximum force during all permutations of stimulus magnitude (mA), duration (us), and temperature:
##
##		ID			Muscle	Stimulus(mA)	Pulse_With(us)	temperature(ºC)	(Max Force)(N)
##		REdC10-M1	M1		10				50				10				result
##		REdC10-M1	M1		20				50				10 				result
##		REdC10-M1	M1		30				50				10 				result
##		.
##		.
##		.
##		REdC10-M1	M1		10				10000			10				result
##		REdC10-M1	M1		20				10000			10				result
##		REdC10-M1	M1		30				10000			10				result
##		.
##		.
##		.
##		REdC10-M1	M1		600				50000			10				result

### right now all the column headers in input files have to be identical i.e. "Stimulus"
##### EDIT THE DIRECTORY ON LINE 58 SO IT JUST USES THE CURRENT DIR

import csv
from decimal import *
import statistics
import re
import os

my_stimulus_print_out=[10,20,30,40,50,60,70,80,90,100,200,300,400,500,600]
my_row_number_list = [33562,66893,100224,133560,166891,200228,233558,266894,300225,333561,366891,400227,433557,466893,500229]

## FIRST STEP: find a list of animals in a folder, read in a list of files for an animal, find all the unique animal names
experiment_list=[]
experiment_dict={}
my_stimulus=()
my_dir='C:/Users/rdelcarlo/Desktop/possibly_useful_for_Sage/sample_rheo'
for filename in os.listdir(my_dir):
	if filename.endswith(".dat.csv"):
		#Filename Example: REdC10-M1-Rheo1000-10C1.dat.csv
		m = re.search('((^\S+?)-(\S+?))-(\D+)(\d+).dat.csv', filename)
		my_experiment=str(m.group(1)) #REdC10-M1
		my_muscle=str(m.group(3)) #M1
		my_stimulus=int(m.group(5)) #1000
		experiment_dict[my_stimulus] = filename
		experiment_list.append(my_experiment)
print("all the unique SnakeID-M# combinations are:")
unique=set(experiment_list)
print(unique)

### SECOND STEP: for each unique SnakeID-M# combination, get the list of files and rank in acending order with respect to both temperature and pulse width
for each_experiment in unique:
	experiment_file_list=[]
	outfile_name=each_experiment + '.max_force_vs_temp_table.csv'
	outfile=open(outfile_name, 'w')
	outfile.write("SnakeID,Muscle,Stimulus(mA),Pulse_With(us),Max_Force(g)\n")
	print("unique experiment" + each_experiment)
	for filename in os.listdir(my_dir):
		if filename.endswith(".dat.csv"):
			if each_experiment in filename:
				print(filename) #REdC11-M1-Rheo10000-35C.dat.csv
				m = re.search('(^\S+?)-(\S+?)-(\D+)(\d+)', filename)
				my_number=int(m.group(4))
				experiment_file_list.append(my_number)
				experiment_dict[my_stimulus] = filename #my_stimulus is defined in line 56, and my_temp() in line 57
	unique_pulse_width=[]
	unique_pulse_width = set(experiment_file_list)
	list_unique_pulse_width=list(unique_pulse_width)
	list_unique_pulse_width.sort()
	print(list_unique_pulse_width)
	for pulse_width in list_unique_pulse_width:
		pulse_width = str(pulse_width)
		muscle_number = each_experiment.split("-")
		snake = muscle_number[0]
		musc_num = muscle_number[1]
		baf = str(each_experiment + "-Rheo" + pulse_width + ".dat.csv")
		if(os.path.isfile(baf)):
			with open(baf) as csvfile:
				print(baf)
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
					print("Sum_stim is greater than zero")
					line_pass = 10000
					for each in stim:
						row_num = row_num + 1
						if each > 2:
							if row_num > line_pass:
								line_pass = row_num
								location = stim.index(each)-1
								section = location -100
								force_list = force[section:location]
								average = (sum(force_list)/100)
								max_force_list = force[(location+1):(location + 1501)]
								length_max_list = len(max_force_list) #not sure why length is important yet
								maximum_force = max(max_force_list)
								table_value = maximum_force - average #table value is the maximum force reported in grams (g), independent of muscle mass
								my_stim_num = my_stimulus_print_out[stimulus_num]
								my_print_list = [snake,musc_num,my_stim_num,pulse_width,table_value]
								for each in my_print_list:
									outfile.write(str(each))
									outfile.write(',')
								print(snake,musc_num,my_stim_num,pulse_width,table_value)
								outfile.write('\n')
								line_pass = line_pass + 1000
								stimulus_num=stimulus_num+1
				elif sum_stim == 0:
					print("sum_stim equals zero")
					list_pos=0
					for each in stim:
						row_num = row_num + 1
						list_len = len(my_row_number_list)
						if list_pos < list_len:
							if row_num == my_row_number_list[list_pos]:
								location = row_num - 1
								section = location -100
								force_list = force[section:location]
								length=len(force_list)
								average = (sum(force_list)/100)
								max_force_list = force[(location+1):(location + 1501)]
								length_max_list = len(max_force_list)
								maximum_force = max(max_force_list)
								table_value = maximum_force - average
								my_stim_num=my_stimulus_print_out[stimulus_num]
								print(snake,musc_num,my_stim_num,pulse_width,table_value)
								my_print_list = [snake,musc_num,my_stim_num,pulse_width,table_value]
								for each in my_print_list:
									outfile.write(str(each))
									outfile.write(',')
								outfile.write('\n')
								list_pos = list_pos + 1
								stimulus_num=stimulus_num+1
