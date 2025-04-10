# tcouchii-tradeoff

This github repository includes all the code and (eventually) all the data required to recreate the analyses and figures in Babish et al ?2025?, which examines the trade-offs in sodium channel, muscular, and whole-animal functioning associated with TTX resistance in the Sierra garter snake, *Thamnophis couchii*. The organization of the files is such that all scripts can be run without any rearranging, and all outputs will end up in the same directory that the current output files are in (overwriting them, so be careful of that if you want to save the originals).

## Myography

This folder contains all the data and scripts associated with the muscle-level analyses presented by Babish et al. This includes C4P, a transient contraction experiment, tetanus, a tetanic contraction experiment, and Rheobase, an excitability experiment. Within folders (data_raw, OutFiles, etc), files are organized into subfolders according to which of these three experiments they fall under. They are organized in this way so that people who are not interested in the raw data but want the outputs, for example, can download just that information more easily.

### data_raw

#### C4P

All the C4P data files are named with the Snake CollPreNo ID, the muscle number (some snakes had multiple muscles dissected out and tested on, and the pulse number (1 through 4). This naming is essential to the functioning of the script 01_C4P_dataprocessing.R and should not be changed. The data has 4 columns: Time (s), Force (g), Position (mm), and Stimulus (?).

#### Tetanus

All the Tetanus data files are named with the Snake CollPreNo ID, the muscle number (some snakes had multiple muscles dissected out and tested on, and which repeat of the protocol it was (1 or 2). This naming is essential to the functioning of the script 04_Tetanus_dataprocessing.R and should not be changed. The data has 4 columns: Time (s), Force (g), Position (mm), and Stimulus (?).

#### Rheobase

#### Snake_data_sheets

These data sheets contain various information about the snakes used in all three of these experiments, and are in this folder because all three analyses draw on this information. The most important information in these files is snake MAMU and muscle masses for the muscles used in the various analyses, the former of which is integral to the analyses and the latter of which is necessary for all data processing (scaling results to muscle mass).

#### IC50

IC50 data, currently unorganized except IC50.csv, which has the final results for all the IC50 analyses and isn't really raw data but is staying there until I have the code to create it and can actually treat it as an outfile.

### data_processed

This folder might end up deleted for all but rheobase; more likely, the force and force-1d files may end up here because they're more data than results and so should be in a data folder, but that's tbd.

### OutFiles

#### C4P

-   Couchii_C4P_Force.csv contains the force output of the muscles, scaled to muscle mass and normalized to baseline

-   Couchii_C4P_Force_Sorted.csv contains the same information as the above file, just sorted and rotated to be better used by the later analyses (except actually I don't think I need it? needs confirmed)

-   Couchii_C4P_Force_1d.csv contains the first derivative of the force data (from every 200th reading)

-   Couchii_C4P_Metrics.csv contains the key summary stats for each transient contraction, including base force, contraction amplitude, contraction duration, and more

-   P1-C4PFiles.csv contains a list of the raw data files that contributed to the above files; it is not directly used for anything but may prove informative or useful

#### Tetanus

-   Couchii_Tetanus_Force.csv contains the force output of the muscles, scaled to muscle mass and normalized to baseline

-   Couchii_Tetanus_Force_Sorted.csv contains the same information as the above file, just sorted and rotated to be better used by the later analyses (except I may not actually use it, need to confirm)

-   Couchii_Tetanus_Force_1d.csv contains the first derivative of the force data (from every 200th reading)

-   Couchii_Tetanus_Metrics.csv contains the key summary stats for each tetanic contraction, including base force, contraction amplitude, contraction duration, and more

-   P1-TetanusFiles.csv contains a list of the raw data files that contributed to the above files; it is not directly used for anything but may prove informative or useful

#### Rheobase

### scripts

-   01_C4P_dataprocessing.R - This file takes all the raw C4P data files (from data_raw/C4P) and processes them to determine amount of force exerted at each time point (normalized and scaled to muscle mass). It then calculates several metrics of interest for the transient contractions, including contraction amplitude, time to maximal contraction, etc., and generates a first derivative trace of the force curve (giving force/second) to determine the minimal and maximal rate of force change. The output files are used by file 02_C4P_analysis.R. Will eventually also remove outliers, still trying to figure out what method I feel comfortable with there.

-   02_C4P_analyses.R - Does correlation analyses (Kendall's tau and Pearson's correlation) between C4P metrics and both MAMU and IC50. Also does basic linear regressions and calculates RMSE and R\^2. Creates .csvs with the correlation values (Couchii_C4P_X_corr_split.csv) and the regression parameters (Couchii_C4P_X_lm.csv) in OutFiles.

-   03_C4P_visualization.R - Creates all C4P-related (publishable) figures. Includes boxplots comparing contraction metrics across pulses, line plots comparing them at an individual level across pulses, trace plots of the contraction force and force derivative broken up by pulses, and the same trace plots with all 4 pulses aggregated.

-   04_Tetanus_dataprocessing.R - This file does basically the same thing as script #1, but for the tetanic contraction protocol data (Tetanus) instead of the transient contraction protocol data (C4P). The output files are used by file 05_Tetanus_analysis.R. Will eventually remove outliers as well/

-   04_Tetanus_analyses.R - Does correlation analyses (Kendall's tau and Pearson's correlation) between tetanus metrics and both MAMU and IC50. Also does basic linear regressions and calculates RMSE and R\^2. Creates .csvs with the correlation values (Couchii_Tetanus_X_corr.csv) and the regression parameters (Couchii_Tetanus_X_lm.csv) in OutFiles.

-   06_Tetanus_visualization.R - Creates all tetanus-related figures, including traces of the force and the first derivative of the force. Ideally it will eventually plot the traces as fitted to the sigmoidal equation from REdC, but that isn't working right now.

-   Absolute nightmare that is the rheobase scripts will be dealt with eventually, after the above and the whole-animal stuff (you can tell I don't want to deal with them right now)

### IC50_analyses

This will ideally eventually be spread out between other folders, but it's all in this folder right now because I'm missing some of the code associated with it and therefore can't do anything with it yet anyway.

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
