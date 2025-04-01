#read in megareport.csv
megareport <- read.csv("possibly_useful_for_Sage/Megareport_to_pw_sorted_PERL/megareport.csv")

#get list of unique pulse values 
pulse_u <- sort(unique(megareport$Pulse_Width_us))
#and list of unique stimulus values, sorted numerically
stim_u <- sort(unique(megareport$Stimulus_mA))

#make list of first two columns (ID and muscle number) with hyphens
 #(will be column headers later)
  #currently also including temp number (meaningless) to make sure other code runs
id_list <- unique(paste(megareport$SnakeID,"_",megareport$Muscle,"_10",sep=""))

#add column so this is actually useful for searching
megareport$IndID <- paste(megareport$SnakeID,"_",megareport$Muscle,"_10",sep="")

#split megareport into list of dataframes by pulse width
dfs <- split( megareport, f = megareport$Pulse_Width_us)

#loop through all pulse widths
for (i in 1:length(pulse_u)) {
  #pull out df with proper pulse width
  pulsedf <- dfs[[i]]
  
  #make output dataframe
  output <- data.frame(matrix(ncol=0,nrow=16))
  #add on Dose column (with 0 for y-intercept)
  output <- cbind(output, Dose = c(0,stim_u))
  
  #initialize vector to assign column names
  colnames <- c("Dose")
  
  #loop through each individual
  for (ind in id_list) {
    #pull out max force from each stimulus value
    forces <- pulsedf[pulsedf$IndID == ind,6]
    
    #some of them don't have all the data for some reason, so ditch those
    if (length(forces) == 15) {
      #find max force
      maxforce <- max(forces)
      #scale all forces according to that max force
      forces_scaled <- forces/maxforce
      
      #add a (0,0) point to the beginning of the vector
      forces_vector <- c(0,forces_scaled)
      
      #add this specimen's data to output dataframe
      output <- cbind(output,ind = forces_vector)
      
      #save colname to fix later (because it's hard to use variables for this)
      colnames <- append(colnames, ind)
      }
  }
  #set column names appropriately
  colnames(output) <- colnames
  
  #for debugging
  print(pulse_u[i])
  
  #write output file (into OutFiles directory, which the next code chunk wants)
  filename = paste("OutFiles/TPA_",pulse_u[i],".csv", sep="")
  write.csv(output, file = filename)
}

#output:
#files TPA_pulse.csv
#first column is unique stimulus values
#subsequent columns are max force of each individual at those different
 #levels of stimulus, scaled to max out of them all (so ranging from 0-1)