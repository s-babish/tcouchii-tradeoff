# tcouchii-tradeoff

This github repository includes all the code and (eventually) all the data required to recreate the analyses and figures in Babish et al (in prep), which examines the trade-offs in sodium channel, muscular, and whole-animal functioning associated with TTX resistance in the Sierra garter snake, *Thamnophis couchii*. The organization of the files is such that all scripts can be run without any rearranging, and all outputs will end up in the same directory that the current output files are in (overwriting them, so be careful of that if you want to save the originals).

## Myography

This folder contains all the data and scripts associated with the muscle-level analyses presented by Babish et al. This includes C4P, a transient contraction experiment, tetanus, a tetanic contraction experiment, and Rheobase, an excitability experiment. Within folders (data_raw, OutFiles, etc), files are organized into subfolders according to which of these three experiments they fall under. They are organized in this way so that people who are not interested in the raw data but want the outputs, for example, can download just that information more easily.

### data_raw

All the raw data for the myography portion of the paper; the only processing before this was converting the files from those used by the proprietary software on our data logger to .csv files.

#### C4P

All the C4P data files are named with the Snake CollPreNo ID, the muscle number (some snakes had multiple muscles dissected out and tested on), and the pulse number (1 through 4). This naming is essential to the functioning of the script 01_C4P_dataprocessing.R and should not be changed. The data has 4 columns: Time (s), Force (g), Position (mm), and Stimulus.

#### Tetanus

All the Tetanus data files are named with the Snake CollPreNo ID, the muscle number (some snakes had multiple muscles dissected out and tested on), and which repeat of the protocol it was (1 or 2). This naming is essential to the functioning of the script 04_Tetanus_dataprocessing.R and should not be changed. The data has 4 columns: Time (s), Force (g), Position (mm), and Stimulus.

#### Rheobase

All the rheobase data files are named with the Snake CollPreNo ID, the muscle number (some snakes had multiple muscles dissected out and tested on), and what the pulse width was (50, 100, 200, 500, 1000, 2000, 5000, 10000, and 50000, in varying permutations for each individual). There is a "1" at the end of each of these numbers for some reason, so that files from pulse width 50 (microseconds) actually end "501.dat.csv"). The data has 4 columns: Time (s), Force (g), Position (mm), and Stimulus.


#### Snake_data_sheets

These data sheets contain various information about the snakes used in all three of these experiments, and are in this folder because all three analyses draw on this information. The most important information in these files is snake MAMU and muscle masses for the muscles used in the various analyses, the former of which is integral to the analyses and the latter of which is necessary for all data processing (scaling results to muscle mass).

#### IC50_analysis

Currently empty, this will contain the data files needed in the calculation of the IC50s for muscle=level TTX resistance once I find them (or I might delete this and stick with IC50.csv in the data_processed folder which is where I've been pulling that data from this whole time anyway).

### data_processed

In theory, this should contain data at an intermediate stage of processing, not ready for analysis or plotting but not raw. This folder may eventually be deleted for all protocols but rheobase; more likely, the force and force-1d files may end up here because they're more data than results and so should be in a data folder. 

### OutFiles

#### C4P

This folder is split into subfolders arranging individuals by species/genotype. Each folder follows the same basic structuring of the file names, just switching out the genotype/species indicator term. 

The contents of each subfolder are:
-   X_C4P_Force.csv contains the force output of the muscles, scaled to muscle mass and normalized to baseline

-   X_C4P_Force_1d.csv contains the first derivative of the force data

-   X_C4P_Metrics.csv contains the key summary stats for each transient contraction, including base force, contraction amplitude, contraction duration, and more

-   X_P1_C4PFiles.csv contains a list of the raw data files that contributed to the above files, along with the snake IDs and muscle masses

Additionally, this folder contains a .csv titled "All_geno_avgs.csv" that contains the mean and standard error for the twitch/transient contractions for each genotype, which can be drawn from to plot. The plotting scripts do not currently draw on it but will once data processing is complete. There are also two .csv files, "Geno_dunn_comparisons.csv" and "Subset_geno_dunn_comparisons.csv" which contain the results for Dunn's post hoc tests for either all the groups or a subset of the groups (WT sirtalis, atratus, hammondii, and LVNV, EPN, and T). 

#### Tetanus

This folder is very similar to the C4P folder previously described, with substructuring to present the results of the tetanic contraction experiments by genotype/species.

The contents of each subfolder are:
-   X_Tetanus_Force.csv contains the force output of the muscles, scaled to muscle mass and normalized to baseline

-   X_Tetanus_Force_1d.csv contains the first derivative of the force data (from every 200th reading)

-   X_Tetanus_Metrics.csv contains the key summary stats for each tetanic contraction, including base force, contraction amplitude, contraction duration, and more

-   X_P2_TetanusFiles.csv contains a list of the raw data files that contributed to the above files; it is not directly used for anything but may prove informative or useful

Additionally, this folder contains a .csv titled "All_geno_avgs.csv" that contains the mean and standard error for the twitch/transient contractions for each genotype, which can be drawn from to plot, and another version "All_geno_avgs_no_outliers.csv" which contains the same data but calculated without observations identified as bad (muscle tearing or twitch rather than transient contractions). The plotting scripts do not currently draw on these files but will once data processing is complete. There are also two .csv files, "Geno_dunn_comparisons.csv" and "Subset_geno_dunn_comparisons.csv" which contain the results for Dunn's post hoc tests for either all the groups or a subset of the groups (WT sirtalis, atratus, hammondii, and LVNV, EPN, and T). 

#### Rheobase

Again substructured to present the results of the rheobase contraction experiments by genotype/species.

The contents of each subfolder are:

- megareport.csv contains the maximum force (not scaled) produced by each muscle for each stimulus level and pulse width.
- Rheobase_#us.csv contain the subsets of the megareport for each pulse width.
- Force_scaled/Rheo#.csv contains the force produced by each muscle, scaled to the maximal force that muscle produced, at each pulse width.
- Sigmoidal/Sigmoidal-rpt.csv contains the full data from the sigmoidal curves fitted to the force-current relationship at each pulse width. The other files in this folder are subsets with individual pulse widths and with (-rpt) or without () further information about the individual the data was collected from.
- Sigmoidal_sub/ contains all the above information but calculated from sigmoidal curves fitted only to the contractions at 0-100mA (although the forces used to make those curves were still scaled to their overall maximal value)


### scripts

-   00_sample_mapping_comparison.R - This script may not last, but for now it exists to map the samples we used in the various myography experiments (I wanted an idea of where they came from) and compiles a list of which samples were used for which experiments (because annoyingly there were very different individuals used for each experiment and there's not huge overlap)).

-   01_C4P_dataprocessing.R - This file takes all the raw C4P data files (from data_raw/C4P) and processes them to determine amount of force exerted at each time point (normalized and scaled to muscle mass). It then calculates several metrics of interest for the transient contractions, including contraction amplitude, time to maximal contraction, etc., and generates a first derivative trace of the force curve (giving force/second) to determine the minimal and maximal rate of force change. The output files are used by file 02_C4P_analysis.R. 

-   02_C4P_analyses.R - Does correlation analyses (Kendall's tau and Pearson's correlation) between C4P metrics and both MAMU and IC50. Also does basic linear regressions and calculates RMSE and R\^2. Finally, performs Dunn's tests on the various C4P metrics between genotypes to see how genotypes vary in their transient contractions; the results from these are stored in "Geno_dunn_comparisons.csv." The subset versions are those that have been performed on datasets with outliers manually removed. 

-   03_C4P_visualization.R - Creates all C4P-related (publishable) figures. Includes boxplots comparing contraction metrics across pulses, line plots comparing them at an individual level across pulses, trace plots of the contraction force and force derivative broken up by pulses, and the same trace plots with all 4 pulses aggregated.

-   04_Tetanus_dataprocessing.R - This file does basically the same thing as script #1, but for the tetanic contraction protocol data (Tetanus) instead of the transient contraction protocol data (C4P). The output files are used by file 05_Tetanus_analysis.R. 

-   05_Tetanus_analyses.R - Does correlation analyses (Kendall's tau and Pearson's correlation) between tetanus metrics and both MAMU and IC50. Also does basic linear regressions and calculates RMSE and R\^2. Finally, performs Dunn's tests on the various tetanus metrics between genotypes to see how genotypes vary in their transient contractions; the results from these are stored in "Geno_dunn_comparisons.csv." 

-   06_Tetanus_visualization.R - Creates all tetanus-related figures, including traces of the force and the first derivative of the force.

-   P1_Rheobase_processing.py takes in the raw rheobase data files and processes them to an intermediate (stored in data_processed)

-   P2_Rheobase_compilation.py takes the output from the previous step and produces the file "megareport.csv" for each genotype/group, which goes into that group's Rheobase OutFiles folder.

-   P3_Rheobase_sorting.R takes the megareport and sorts the data into outfiles by pulse length (and also scales the force data for each individual). Outputs for this are stored in the Force_scaled subfolder of the OutFiles for each group, and it may be moved to data_processed eventually. Those outputs have the format "Rheo_pulse width," with pulse widths 10, 50, 100, 200, 500, 1000, 10000, 50000, 100000, and 500000, with some occasionally missing if that data was not recorded on that muscle.

-   X1_Rheobase_sigmoidal.R - This fits a sigmoidal curve to the scaled force data obtained from script P3. Changing which lines are commented and some folder names switches this between working with the full range of currents used in the rheobase protocols and working only with those below 200 mA, which is the value identified by Bobby, Norm, and myself as where the results cease following the expected sigmoidal curve. Results from the former are saved in folders within Outfiles and their respective genotypes labeled "Sigmoidal," while those calculated from a subset of the current values are stored in the "Sigmoidal_sub" folders. This script also plots the data and the fitted curve, and at some point I need to subdivide those plots further and use them to identify problematic observations or muscles.

-   X2_Rheobase_analyses.R - This fits the equation that produces the rheobase strength-duration curve and then performs statistics to compare the rheobase and chronaxie parameters across genotypes. In the process, it plots the raw data and the fitted curves so that you can ensure the curves are being fitted properly. It also performs stats on the various sigmoidal curve parameters, include x0, dx, ec10, ec50, and ec90. Those stats are still done across pulse lengths instead of genotypes; that needs to be updated. 

-   X3_Rheobase_visualization.R - Makes a lot of plots, I need to go through and delete the ones I know I won't want to use. Current highlight is the plotting of the average fitted strength-duration curves for each genotype, with potential to add a ribbon with standard error. Also made the rheobase plots in the supplemental slides of the JMIH talk.

-   The "unused" folder that contains things I got from Bobby that I don't think I need but don't want to delete just yet

### Plots
This folder was made to compile all the plots from subfolders within OutFiles, in a somewhat more organized manner, and with a focus on plots that are more likely to end up in the eventual manuscript. They are still sorted into subfolders by protocol.

## tradeoff_modeling

This folder contains all the data and code associated with the whole-animal trade-off analyses performed by Babish et al. Miscellaneous research design scripts I plan to pull from and eventually remove are just chucked in the top level of the folder at the moment.

### data_raw

This file contains the data files the rest of the project draws from. couchiionlydata.csv contains old data collected by EDBIII and his lab, and some old data collected by CRF and his lab. summer24couchii.csv contains the additional data collected by SDB and CRF in 2024. old_MAMUs.csv is used to populate MAMU data for some of the individuals that are missing MAMUs in the main datafile(s). couchiiwithoutresistance.csv and fastcouchiineonates.csv were reference files that will probably be deleted, and multispeciesdata.csv was used for a different version of this modeling and may be removed (unless we decide to include that information in the manuscript for comparison).

### OutFiles

This file will have data with information about models, maybe? For now it exists as a just-in-case.

### scripts

-   01_datamanipulation.R - This script merges the raw data frames, scales the data in preparation for modeling, and does some data exploration that was used to determine the models used in script 2. I may or may not do some imputing as well, I want to ask some people with more experience about it first.

-   02_modeling.R - This script runs GLMs, GLMMs (eventually, once I fill in locality information and standardize it), and path models on the whole-animal trade-off data.

-   03_visualization.R - This script will eventually make whatever plots I want for this section, I haven't decided what those look like yet so for now it's empty.

## Ephys

This folder contains all data, and eventually scripts and figures, related to the patch-clamp electrophysiology work that is assessing what effects TTX resistance-conferring mutations have on Nav1.4 function.

### Bobby_masterdocs_results 
This folder contains the masterdocs of patch-clamp electrophysiology results that were used in del Carlo et al 2024. These are not meant to be edited, and they are here mostly in case I mess something up while appending my own electrophysiology results.

### Appended_masterdocs_results
This folder contains masterdocs of patch-clamping data/results with the data I have collected/derived appended on to the data from del Carlo et al 2024. These will be used for analysis later.

### data_raw
The usual drill, this folder contains raw data, sorted by genotype and cell identity. Each *.abf file is named with the date, genotype, cell number, and protocol. Not all cells have data for all protocols, as seals often ruptured before every protocol could be run.

## Presentations

This folder contains any and all powerpoints and other presentations made for this project. This will likely be removed when closer to publication but is useful for now to help share materials, previous outlines, etc. between coauthors. List of presentations below:

- Babish_JMIH 2025 is slides from the 12 minute talk given at JMIH 2025 in ST. Paul, MN presenting preliminary results
- Figure powerpoint that contains the set of plots we are more likely to use in the final publication, with explanations in the descriptions. All plots here can also be found in the "Plots" tab under Myography (the reverse is not true).
