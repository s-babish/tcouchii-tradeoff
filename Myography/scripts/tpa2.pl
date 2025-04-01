#!/usr/bin/perl

# This script will read in the "megareport.csv". 
# We will create unique IDs that include ID_Muscle_Temp, we want the result to be a list of MaxForce values for different Stimulus values
# We will have multiple outfiles (which we will send to a directory) that lump measurements based on a single Pulse value

# on the command line IN THE DIRECTORY WHERE YOUR MEGAREPORT.CSV FILE IS! make sure you have a directory named "OutFiles"
# type on command line :perl tpa.pl megareport.csv

use strict; use warnings;

my $in     = $ARGV[0];			# This will be the "megareport.csv" file that we want to read in, specified on the command line
my $files  = "OutFiles";		# This is the name of a premade directory where we want to store all our outfiles

open (IN, "< $in") or die "Couldn't open IN file\n";	# Open the file on the command line or stop running the script if you can't open it!
<IN>;												#strips off column names/first line of IN file aka megareport.csv

my %pul;							#declaring a bunch of hashes
my %uniqs_l;
my %uniqs_s;
my %stim_lists;

my @lin;

my $uniq_l;
my $uniq_s;
my $eff;
my $us = "_";								# This little buddy is here to save the day a little later on...because numbers in regexes are real dumb

while (<IN>) {								# Here we loop through our IN file line by line
	chomp;									# Always use 'chomp' to get rid of line endings
	@lin				 = split(",", $_);  # for each line, split its contents by a comma
	$pul{$lin[3]}		 = 1;				# we are making a hash of pulse values to get unique list of values
	$stim_lists{$lin[2]} = 1;				# hash of unique stimulus values 
	$uniq_s				 = join("_", $lin[0], $lin[1], $lin[4], $lin[3], $us);  # here we are making a unique identifier composed of ID, Muscle, temp, and pulse...with two extra rogue underscore (FTW!) to keep number matching in line!
	$uniq_l				 = join("\ ", $uniq_s, $lin[2]);		# Second unique identifier that adds the stimulus value to the end of above $unique variable/scalar
	$lin[5]				 =~ m/(\S+)\s/;			# Turns out that our MaxForce value held on to a weird whitespace character when we split by commas, so we are separating the value from that here
	$eff				 = $1;					# This is our MaxForce value w/o weird whitespace value....I was frustrated, hence the uninformative scalar name...


	$uniqs_l{$uniq_l} = $eff; 		# Make a hash of the second unique id we made as the keys and the MaxForce value measured as the values
	$uniqs_s{$uniq_s} = 1;			# Unique hash of the first unique id we made :]
}

close(IN);							# Closing out our IN file



my %id_pul_force;					# Declaring a new hash

my @pulses = keys %pul;					# Extracting the unique pulse values from the hash we made earlier
my @stim_list = sort {$a <=> $b} keys %stim_lists; # Pulling the stimulus values and sorting them numerically (the {$a <=> $b} part, since perl defaults to sorting alphabetically

my $stimuli = join(",", @stim_list);		# Turning our sorted stimulus values into a scalar with values separated by commas for the header of our outfiles

foreach my $short (keys %uniqs_s) {			# Looping through the keys of the %uniqs_s hash that has MaxForce values
	my %stim;
	
	my @forces;
	my @longs;
	
	foreach my $long (keys %uniqs_l) {			# this says to (within the foreach loop above) also loop through the %uniqs_l hash
		@longs = split("\ ", $long);			# We are splitting the keys of %uniqs_l hash by the space (which has the stimulus value as the second part)
		if ($short =~ m/($longs[0])/) {			# Now we are matching the first part of the $long keys (everything but the stimulus) with the %uniqs_s keys
			$stim{$longs[1]} = $uniqs_l{$long}; # essentially we are making a hash that has each stimulus value (for the unique id/musc/temp/pulse combo) as keys and MaxForce as the values
		}
	}
	foreach my $stimulus (sort {$a <=> $b} keys %stim) { # now we sort the hash we just made so the MaxForce values will match up with the ascending stimulus values in our outfiles
		push(@forces, $stim{$stimulus});				 # we add MaxForce values in the proper order to an @array (here @forces)
	}
	
	my $force = join(",", @forces);						# This turns @forces into a $scalar so the values are separated by commas
	
	$id_pul_force{$short} = $force;						# we now have a hash where the keys are the id_musc_temp_pulse and the values are the list of observed MaxForces (in $scalar form, comma separated)
}
	
foreach my $pulse (@pulses) {					# This is where we separate everything out by pulse value so we can write out MaxForce measurements to different files based on pulse values
	open (OUT, "> $files//TPA_$pulse\.csv") or die "Couldn't open OUT for $pulse\n"; # We loop through the list of pulse values and make new files for each that are sent to our OutFiles/ directory
	print OUT "ID\,$stimuli\n";					# Printing the header line

	foreach my $ids (keys %id_pul_force) {			# Simply match the pulse in the unique id in %id_pul_force hash from lin 73 and print out our unique id (now '-' separated) and our list of MaxForce values!
		my @id = split('_', $ids);
		if ($id[3] == $pulse) {
			my $ID = join('-', $id[0], $id[1], $id[2]);
			print OUT "$ID,$id_pul_force{$ids}\n";		# Neat!
		}
	}
	
	close(OUT);										# always close out! Otherwise we will get some weird stuff happening with the outfiles...
}
