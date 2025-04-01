							##!/usr/bin/env python
### Created by Julie Allen Ph.D. and modified by Robert del Carlo (in prep.)
### print out one table per animal with with the maximum force values from 15 locations in the file
### these values will be plotted against known stimulus magnitudes (10(d10)100, 100(d100)600)
### this will be done for six files recorded from each individual; each at a different temperature
### The temperatures are: 10, 15, 20, 25, 30, and 35 (ºC) - they were performed in random order
##
##	filename = REdC10-M1-Rheo10000-35C.dat.csv
##		animal ID = REdC10-M1   ## now include the muscle name
##		muscle = M1
##		protocol = Rheo10000
##		pulse witdth = 10000
##		Temperature = 35(ºC)
##
##		first check that the stimulus column does not = 0
##			if It does then use set row numbers
##			if NOT then read through stimulus column for greater accuracy
##				find the first instance where the stimulus is greater than 2
##				then from 1 less than  position take the average of the PRIOR 100 cells of FORCE column = A
##				then starting at cell where stimulus is > 2 take 1500 cells forward and find the maximum Force in that mix = B
##				Report B - A (maximum force) into the table
##
##  Report one table per animal with:
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
##### EDIT THE DIRECTORY ON LINE 60 SO IT JUST USES THE CURRENT DIR

import csv
from decimal import *
import statistics
import re #regex
import os #operating system - allows it to read files in the directory and interact with the computer (essentially like package in R)

my_stimulus_print_out=[10,20,30,40,50,60,70,80,90,100,200,300,400,500,600]
my_row_number_list = [33562,66893,100224,133560,166891,200228,233558,266894,300225,333561,366891,400227,433557,466893,500229]

## find all the animals in a folder; then for each animal, go through the process to create a table

## FIRST STEP
## read in a list of files for an animal find all the unique animal names
experiment_list=[] #brackets means list
experiment_dict={} #curly brackets means dictionary, a definition for an object (key value pairs)
my_stimulus=() #creates an empty variable

my_dir='C:/Users/rdelcarlo/Desktop/possibly_useful_for_Sage/sample_rheo'
for filename in os.listdir(my_dir):
	#create a list of experiments in the folder
	if filename.endswith(".dat.csv"):
		#print(filename)
		#m = re.search('((^\S+?)-(\S+?))-(\D+)(\d+)', filename)
		#This is the 'regular experession' which parses the filenames by hyphens (SnakeID-M#-Rheo1001.dat.csv)
		#It searches the filenames in my_dir string for... (concat)(enated)
		#SnakeID(my_experiment): (matching ^any character \S except nonwhitespace + repeating one or more times ? the first instance)
		#This ? restricts the "greediness" of the Regular Expression
		#MuscleNumber(my_muscle): (matching \S nonwhitespace characters + repeating one or more times ? optionally)
		#StimulusDuration(my_stimulus): (\D Any nondecimal digit + repeating one or more times)(\d any decimal digit + repeating one or more times)
		m = re.search('((^\S+?)-(\S+?))-(\D+)(\d+)1.dat.csv', filename)
# 		# This adds another component to the filename which consists of repeating decimal digits concatenated with the alphanumeric string C1 as in REdC10-M1-Rheo10000-25C1.dat.csvfile
		my_experiment=str(m.group(1))
# 		#m.group(1)=((^\S+?)-(\S+?)) the first group refers to the SnakeID-M# Combination
		my_muscle=str(m.group(3))
# 		#m.group(3)=(\S+?) the third group refers to the M# string
		my_stimulus=int(m.group(5))
# 		#The fifth group is the intergers in the concatenated combination of Rheo###+
#		#my_stimulus=int((my_indicated - 1)/10)
# 		#This takes the filename artifact ##+1 and converts it to the stimulus pulse width in microseconds e.g. Rheo1001 --> 100 microseconds
# 		#Comment out "my_stimulus" and rename "my_indicated" as "my_stimulus" as filenaming artifact not present in these files
#		my_temp=int(m.group(6))
#		print(filename,my_experiment,my_muscle,my_stimulus,my_temp)
		experiment_dict[my_stimulus] = filename
		#The object experiment_dict{} described above is populated by the newly calculated stimulus pulse width which is renamed here in the filename
		experiment_list.append(my_experiment)
		#The object experiment_list is populated with the SnakeID-M#
print("all the unique SnakeID-M# combinations are:")
unique=set(experiment_list)
# #defines "unique" as a set consisting of all items which meet the requirements of experiment_list and therefore should be unique.
# #consider changing "experiment_list" to "experiment_list" to avoid the SnakeID-M# vs. animal redundancy
print(unique)
# #this will appear as following in the teminal "all the unique animals are {'REdC11-M1', 'Wills06-M1'}"
# #this will appear as following in the teminal "all the unique SnakeID-M# combinations are {'REdC11-M1', 'Wills06-M1'}"
#
# ########  OPEN AN OUTPUT FILE FOR EACH ANIMAL
# #outfile_name=each_experiment +'.table.txt'
# #outfile=open('workfile', 'w')
# ### PRINT THIS TO THE OUTPUT
# #outfile.write("stimulus,location,length,average,length_max,max_force,table")
# #print("stimulus,location,length,average,length_max,max_force,table")
#
# #animals_overall_dict={}
#
# ### SECOND STEP -- for each unique SnakeID-M# combination, get the list of files and rank in acending order.
## for each unique SnakeID-M# combination, get all the files, calculate the stimulus for each animal and order the files by the stimulus
for each_experiment in unique:
	#Create an outfile and write the first line to the outfile (column headers)
	#I believe all the for loops here are creating lists which will then be printed as a list set
	experiment_file_list=[]
	#Defines an object experiment_file_list that will be populated by the results of the for loop
	outfile_name=each_experiment + '.max_force_table.txt'
	#defines the naming convention of the output file
	outfile=open(outfile_name, 'w')
	#the 'w' creates a text file for writing with the stream positioned at the beginning of the file
	#outfile.write("ID\tMuscle\tStimulus (mA)\tPulse_With\t(B-A)\n")
	#\t breaks the column headers into different columns
	outfile.write("SnakeID\tMuscle\tStimulus(mA)\tPulse_With(us)\tMax_Force(N)\n")
	print("unique experiment " + each_experiment)
	#This ^for loop creates a text file for every animal, printing into the terminal "Unique animalSnakeID-M#" followed by the filenames which matched.
	#I believe this vfor loop prints the second column of the output file, populated with the M# written for example M1
	for filename in os.listdir(my_dir):
	#Print the filename for each file in my_dir if the filename ends with .dat.csv and if it belongs to the set each_experiment
	#I believe this ^for loop prints only to the terminal
		if filename.endswith(".dat.csv"):
			if each_experiment in filename:
				#REdC11-M1-Rheo100001.dat.csv
				#REdC11-M1-Rheo10000-35C.dat.csv
				print(filename)
				#this prints into the terminal the filenames of all the files who share the same SnakeID-M# combination
				#m = re.search('(^\S+?)-(\S+?)-(\D+)(\d+)', filename)
				#MTJH402-M1-Rheo10000-35C1.dat.csv
				m = re.search('(^\S+?)-(\S+?)-(\D+)(\d+)', filename)
				my_number=int(m.group(4))
				#Each group is numbered in order from the beginning of the regex, that is, group(4) here is (\d+) which refers to the 1001 in Rheo1001
				#Note that this regex does not follow the same grouping as in line 65 as it is missing the (grouping SnakeID-M#)
				#Printed into the terminal will be "[501, 1001, 2001, 5001, 10001, 100001, 500001]"
				experiment_file_list.append(my_number)
				#Appends the list defined in line 110 with the number pulled from the filename via the my_number regex component
				experiment_dict[my_stimulus] = filename
				#experiment_dict{} was created in line 58 left undefined; my_stimulus in line 59, left undefined, but here, both are said to be in the filename
	unique_pulse_width=[]
	unique_pulse_width = set(experiment_file_list)
	list_unique_pulse_width=list(unique_pulse_width)
	list_unique_pulse_width.sort()
	print(list_unique_pulse_width)
	for pulse_width in list_unique_pulse_width:
		#CRF3066-M3-Rheo1000-10C1.dat.csv
		pulse_width = str(pulse_width)
		baf = str(each_experiment + "-Rheo" + pulse_width + ".dat.csv")
		print(baf)
		if(os.path.isfile(baf)):
			with open(baf) as csvfile:
				print(baf)
				reader=csv.DictReader(csvfile) #define the reader as dictreader to use csv's
				## first check to see that the stimulus line does not sum to zero
				sum_stim = Decimal(0) #define the object sum_stim as the decimal zero
				stim=[] #create an list named stim but leave undefined
				force=[] #create an list named force but leave undefined
				#print(sum_stim=sum(Decimal(row['Stimulus']) for row in reader)
				for row in reader: #see columns named Stimulus and Force (g) and for each row, see that t
					sum_stim += Decimal(row['Stimulus']) #+= iterates over every row in 'Stimulus' column and adds to sum_stim which was said zero; if it's still zero at the end?
					stim.append(Decimal(row['Stimulus'])) #add to the object stim all the rows?
					force.append(Decimal(row['Force (g)'])) #add to the object force all the rows?
				#print(sum_stim)
				row_num=0 #start at the first row, absent headers?
				stimulus_num=0 #define stimulus number as zero
				if sum_stim > 0: #if after all the rows of +=, sum_stim is still zero, then print sum_stim is greater than 0
					print("Sum_stim is greater than zero")
					#print (sum_stim, " ","is greater than 0")
					line_pass = 10000 ## this defines the object "line_pass" as 10000
					##  make sure sum_stim is > 0 then do the following if the entire column is nonzero,
					for each in stim: #working row by row, look for any instance of a value >2, then jump
						row_num = row_num + 1
						if each > 2:
							#####
							##### after adjusting tabs, we appear to be stuck in an infinite loop here - could it be that print statement?
							#####
							if row_num > line_pass: #if the current row number is over 10000 higher, and then
								line_pass = row_num #when the current row number reaches the line pass
								#print(row_num,line_pass) #print that coordinate
								location = stim.index(each)-1 #this is the new frame of reference
								section = location -100 #this defines the earliest point in the set which comprises A (the average )
								force_list = force[section:location] #defines force_list as values in force column in the rows from section to location (100 cells)
								average = (sum(force_list)/100) #this defines A
								max_force_list = force[(location+1):(location + 1501)] #this defines the region where the maximum force will be read
								length_max_list = len(max_force_list) #not sure why length is important yet
								maximum_force = max(max_force_list) #searches for B, the maximum force reached in the region defined as max_force_list
								table_value = maximum_force - average #defines B-A to achieve the maximum force value (in g, in N?)
								#table_value = (maximum_force - average)*0.00980665 #converts grams of force to Newtons of force
								#how to integrate with skeletalmusclemasses.csv to then table_value = ((maximum_force - average)*0.00980665)/(SnakeID-M#_mass)???
								print(stimulus_num)
								my_stim_num = my_stimulus_print_out[stimulus_num] #stimulus_num is defined =0 in line 187
								#print(my_experiment,my_muscle,stimnum,my_stimulus,temp,table_value)
								#this prints the five columns of the outfile, SnakeID-M#, M#, Stimulus Magnitude (mA), Stimulus Pulse Width (microseconds), Max_Force (g)
								#print(my_experiment,my_muscle,stimnum,my_stimulus,my_temp,table_value)
								##this prints the 6 columns of the outfile, SnakeID-M#, M#, StimMag (mA), StimDur (us), Temperature(ºC) Max_Force (g hopefully N)
								#my_print_list = [my_experiment,my_muscle,stimnum,my_stimulus,table_value]
								my_print_list = [each_experiment,my_muscle,my_stim_num,pulse_width,table_value] #should this have \t?
								for each in my_print_list:
									outfile.write(str(each))
									outfile.write(',') # \t or ",", in both cases, the output is printing all in one row rather than over columns
								print(each_experiment,my_muscle,my_stim_num,pulse_width,table_value,row_num)
								#print('\n')
								#outfile.write(my_experiment,my_muscle,stimnum,my_stimulus,table_value)
								#print(each,location,length,average,length_max_list,maximum_force,table_value)
								outfile.write('\n')
								line_pass = line_pass + 1000
								stimulus_num=stimulus_num+1
				elif sum_stim == 0:
					print("sum_stim equals zero")
					list_pos=0
					for each in stim:
						row_num = row_num + 1
						#print(row_num,list_pos,my_row_number_list[list_pos])
						list_len = len(my_row_number_list)
						if list_pos < list_len:
							if row_num == my_row_number_list[list_pos]:
								location = row_num - 1
								#location = stim.index(each)-1
								section = location -100
								force_list = force[section:location]
								length=len(force_list)
								average = (sum(force_list)/100)
								max_force_list = force[(location+1):(location + 1501)]
								length_max_list = len(max_force_list)

								#### MAY BE SPACING ISSUE 'EMPTY LIST' MAX FORCE LIST DOESN'T EXIST
								#As the output stands now, the values are not printing to columns but successive rows

								maximum_force = max(max_force_list)
								table_value = maximum_force - average
								my_stim_num=my_stimulus_print_out[stimulus_num]
								print(each_experiment,my_muscle,my_stim_num,pulse_width,temp,table_value)
								#print(my_experiment,my_muscle,stimnum,my_stimulus,my_temp,table_value)
								my_print_list = [each_experiment,my_muscle,my_stim_num,pulse_width,temp,table_value]
								#my_print_list = [my_experiment,my_muscle,stimnum,my_stimulus,my_temp,table_value]
								for each in my_print_list:
									outfile.write(str(each))
									outfile.write(',')
								outfile.write('\n')
								#outfile.write(my_experiment,my_muscle,stimnum,my_stimulus,table_value)
								#outfile.write(my_experiment,my_muscle,stimnum,my_stimulus,my_temp,table_value)
								#print(each,location,length,average,length_max_list,maximum_force,table_value)
								#print(each,location,length,average,length_max_list,maximum_force,my_temp,table_value)
								list_pos = list_pos + 1
								stimulus_num=stimulus_num+1
