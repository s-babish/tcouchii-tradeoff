# tcouchii-tradeoff

This github repository includes all the code and (eventually) all the data required to recreate the analyses and figures in Babish et al ?2025?, which examines the trade-offs in sodium channel, muscular, and whole-animal functioning associated with TTX resistance in the Sierra garter snake, *Thamnophis couchii*. The organization of the files is such that all scripts can be run without any rearranging, and all outputs will end up in the same directory that the current output files are in (overwriting them, so be careful of that if you want to save the originals).

## Myography

This folder contains all the data and scripts associated with the muscle-level analyses presented by Babish et al. This includes C4P, a transient contraction experiment, tetanus, a tetanic contraction experiment, and Rheobase, an excitability experiment. Within folders (data_raw, OutFiles, etc), files are organized into subfolders according to which of these three experiments they fall under. They are organized in this way so that people who are not interested in the raw data but want the outputs, for example, can download just that information more easily.

### data_raw

#### C4P

All the C4P data are named with the Snake CollPreNo ID, the muscle number (some snakes had multiple muscles dissected out and tested on, and the pulse number (1 through 4). This naming is essential to the functioning of the script 01_C4P_dataprocessing.R and should not be changed. The data has 4 columns: Time (s), Force (g), Position (mm), and Stimulus (?).

#### Tetanus

#### Rheobase

#### Snake_data_sheets

These data sheets contain various information about the snakes used in all three of these experiments, and are in this folder because all three analyses draw on this information. The most important information in these files is snake MAMU and muscle masses for the muscles used in the various analyses, the former of which is integral to the analyses and the latter of which is necessary for all data processing (scaling results to muscle mass).

### data_processed

### OutFiles

#### C4P

-   Couchii_C4P_Force.csv contains the force output of the muscles, scaled to muscle mass and normalized to baseline

-   Couchii_C4P_Force_Sorted.csv contains the same information as the above file, just sorted and rotated to be better used by the later analyses

-   Couchii_C4P_Force_1d.csv contains the first derivative of the force data (from every 200th reading)

-   Couchii_C4P_Metrics.csv contains the key summary stats for each transient contraction, including base force, contraction amplitude, contraction duration, and more

-   P1-C4PFiles.csv contains a list of the raw data files that contributed to the above files; it is not directly used for anything but may prove informative or useful

#### Tetanus

#### Rheobase

### scripts

-   01_C4P_dataprocessing.R - This file takes all the raw C4P data files (from data_raw/C4P) and processes them to determine amount of force exerted at each time point (normalized and scaled to muscle mass). It then calculates several metrics of interest for the transient contractions, including contraction amplitude, time to maximal contraction, etc., and generates a first derivative trace of the force curve (giving force/second) to determine the minimal and maximal rate of force change. The output files are used by file 02_C4P_analysis.R.

### IC50_analyses

This will ideally eventually be spread out between other folders, but it's all in this folder right now because I'm missing some of the code associated with it and therefore can't do anything with it yet anyway.

## tradeoff_modeling

This folder contains all the data and code associated with the whole-animal trade-off analyses performed by Babish et al. (it's unorganized right now but I'm going to fix it I promise, one thing at a time)
