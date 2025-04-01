##!/usr/bin/env python
### Created by Julie Allen Ph.D. and Robert Eugene del Carlo
### Compile all *.rheobase.csv report files from rheobase_analysis.py
##  Input File format:
##
##		ID			Muscle	Stimulus(mA)	Pulse_With(us)	temperature(ºC)	(Max Force)(N)
##		REdC10-M1	M1		10				50				10				result
##		REdC10-M1	M1		20				50				10 				result
##		REdC10-M1	M1		30				50				10 				result
##		.			.		.				.				.				.
##		REdC10-M1	M1		600				50000			10				result
### Compile all the *.rheobase.csv files into a concatenated megareport, then sort by pulsewidth an create subsetted outfiles
### 	filename = p4C4PwithTTX-ContrAmpl.N.g..csv - but we need to make an R code specific for this protocol, like TPA-50us.csv
###		Dose	BRN_M3_10C		BRN_M3_15		SnakeID_M#_T*C
### 1	10		0.12806579		0.1359381			...
### 2	20		0.129739542		0.1360012			...
### 3	30		0.123182693		0.1359005			...
### 4	40		0.120206392		0.1429580			...
###	...	...		...									...
###	...	...		...									...
###	...	...		...									...
### 15	600		1.202063950		1.3987429			...
### Then sort by pulse width and make outfiles subsetted by pulse width (i.e. Rheobase50.csv ... Rheobase50000.csv)
###				1			2			3			4			5	6	7	8	9	10	11	12	13	14	15
###Dose(mA)		10			20			30			40			50	60	70	80	90	100	200	300 400 500	600
###BRN_M3_10C	0.12806579	0.129739542	0.123182693	0.120206392	........................................1.202063950
###................................................................................................................
###................................................................................................................
###BRN_M3_35C	0.1359381	0.1360012	0.1359005	0.1429580	........................................1.398742901
###................................................................................................................
### This output format will need to be transposed for use in pipeline linking to R codes for sigmoidal analysis (best by Datamash transpose)
##### EDIT THE DIRECTORY ON LINE 50 SO IT JUST USES THE CURRENT DIRECTORY WHERE DATA FILES AND CODE ARE STORED

import csv
from decimal import *
import statistics
import re
import os
import pandas as pd
import numpy as np

my_stimulus_print_out=[10,20,30,40,50,60,70,80,90,100,200,300,400,500,600]
my_row_number_list = [33562,66893,100224,133560,166891,200228,233558,266894,300225,333561,366891,400227,433557,466893,500229]
my_out_row_number_list = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

experiment_list=[]
experiment_dict={}
my_stimulus=()
my_temp=()
my_dir='C:/Users/sdbab/Downloads/Tcouchii_side_hustle/Tcouchii_side_hustle/Couchii_Rheobase'

skip = '\n'
output_list=[]
for output in os.listdir(my_dir):
	if output.endswith(".rheobase.csv"):
		#Filename Example: CRF3072-M2.max_force_vs_temp_table
		output_list.append(output)
print("All the output files are: ")
OLS=set(output_list)
print(OLS,skip)
snake_list=[]
muscle_list=[]
stim_list=[]
PW_list=[]
for each_report in OLS:
	if(os.path.isfile(each_report)):
		report = pd.read_csv(each_report, header=0, index_col=False)
		snake = set(report['SnakeID'])
		snake_list.append(snake)
		musc = set(report['Muscle'])
		musclist = list(musc)
		muscle_list.extend(musclist)
		stim = set(report['Stimulus(mA)'])
		stimlist = list(stim)
		stim_list.extend(stimlist)
		pw = set(report['Pulse_Width(us)'])
		PW_list.extend(pw)
snake_list = list(filter(None, snake_list))
snake_list.sort()
print("All snakes are: ")
print(snake_list,skip)
muscle_list.sort()
ml = np.array(muscle_list)
uml = np.unique(ml)
uml.tolist()
uml = list(filter(None, uml))
print("Muscle numbers used in this set are: ")
print(uml,skip)
stim_list.sort()
sl = np.array(stim_list)
usl = np.unique(sl)
usl.sort()
usl.tolist()
usl = list(filter(None, usl))
print("All stimuli applied in these experiments are (in mA): ")
print(usl,skip)
PW_list.sort()
pwl = np.array(PW_list)
upwl = np.unique(pwl)
upwl.sort()
upwl.tolist()
upwl = list(filter(None, upwl))
print("All pulse widths used in these experiments are (in microseconds): ")
print(upwl,skip)

pw_sorted_file_list=[]
SIL=[]
analysis_order=[]
#upwl.tolist()
for PW in upwl:
	PW=str(PW)
	pwfile_name = 'Rheobase_' + PW + 'us.csv'
	pwfile=open(pwfile_name, 'w')
	pw_sorted_file_list.append(pwfile_name)
	my_out_row_number_list=str(my_out_row_number_list)
	my_stimulus_print_out=str(my_stimulus_print_out)
	header1 = str(my_out_row_number_list)
	header1 = header1.strip('[]')
	header2 = "SnakeID_Muscle," + my_stimulus_print_out.strip('[]')
	headers = header1 + "\n" + header2
	pwfile.write(headers)
	for each_report in OLS:
		read = pd.read_csv(each_report, header=0, index_col=False)
		snake = read['SnakeID']
		SIL.extend(snake)
USIL=set(SIL)
USIL=list(USIL)
USIL.sort()
print("The Unique snake list is: ")
print(USIL,skip)
print("The total number of snakes in this analysis is: ")
print(len(USIL),skip)

for i in USIL:
	for i in OLS:
		if (i not in analysis_order):
			analysis_order.append(i)
		elif (i in analysis_order):
			continue

print("The list of reports to compile is: ")
print(OLS,skip)
print("The total number of reports to compile is: ")
print(len(OLS),skip)
analysis_order.sort()
print("The compilation order is: ")
print(analysis_order,skip)
print("The total number of files to compile is: ")
print(len(analysis_order),skip)

megareport=open('megareport.csv', 'w')
megareport.write("SnakeID,Muscle,Stimulus_mA,Pulse_Width_us,Max_Force_g\n")
for j in analysis_order:
	report = pd.read_csv(j, header=0, index_col=False)
	report.to_csv(megareport, header=0, index=False, lineterminator='\n')

# print("All reports have been compiled into the megareport. The data will be sorted by pulse width (us) and reported in the following files: ")
# print(pw_sorted_file_list,skip)

# print("The head of the megareport is as follows: ")
# megareport = pd.read_csv("megareport.csv", names=['SnakeID','Muscle','Stimulus_mA','Pulse_Width_us','Max_Force_g'], header=0, index_col=False)
# print(megareport.head(),skip)
# print("For each pulse width (us), rheobase_analysis_compiler will read megareport and sort its contents into files by the following pulse widths: ")
# print(upwl)

# stim_sorted_Fmax=[]
# for k in pw_sorted_file_list:
# 	if(os.path.isfile(k)):
# 		n = re.search('Rheobase_(\\d+)us.csv',k)
# 		j = int(n.group(1))
# 		megareport = pd.read_csv("megareport.csv", names=['SnakeID','Muscle','Stimulus_mA','Pulse_Width_us','Max_Force_g'], header=0, index_col=False)
# 		megareport[megareport.Pulse_Width_us == j]
# 		megareport['SnakeID-Muscle'] = megareport['SnakeID'].str.cat(megareport['Muscle'],sep="-")
# 		megareport.Temperature_C = megareport.Temperature_C.astype(str)
# 		megareport['SnakeID-Muscle-Temperature'] = megareport['SnakeID-Muscle'].str.cat(megareport['Temperature_C'],sep="-")
# 		experiments = np.unique(megareport['SnakeID-Muscle-Temperature'])
# 		experiments.tolist()
# 		for each in experiments:
# 			for farkel in megareport['SnakeID-Muscle-Temperature']:
# 				if each == farkel:
# 					mini = megareport[['SnakeID-Muscle-Temperature','Stimulus_mA','Max_Force_g']]
# 					stim_sorted_Fmax.append(mini)
# print(stim_sorted_Fmax)
