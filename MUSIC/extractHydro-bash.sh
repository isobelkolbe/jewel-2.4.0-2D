#!/bin/bash
# file name: extractHydro.sh
# Written by Isobel Kolbe, 2023
# This script will use a slightly modified version of Chun Shen's `hydro_analysis`
# script to extract temperature and velocity profiles for use with jewel-2.4.0-2D.

# This script needs pandas
# Sometimes the MUSIC directories are structured differently.  If the extractor cannot
# find the data file, try uncommenting line 33.

# Location of the hydro analysis code that will do the extraction
# musicRepoDirectory=$PWD/../hydro_analysis
musicRepoDirectory=~/Documents/Research/JEWEL/myMUSIC/hydro_analysis


shopt -s extglob
# Loop over all the run directories
dirList=(`ls -d */`)
i=0


for dir in $dirList
do 
    echo $dir
    
    dataFile=$dir/evolution_all_xyeta.dat
    if test -f "$dataFile"; then
        cd $dir 
        echo "In " $PWD 
         
        
        mv NcollList* NcollList.dat
        # mv evolution_all_xyeta.dat results/   #If the extractor cannot find the datafile, try uncommenting this line.
        # Clean the directory
        rm -v !("evolution_all_xyeta.dat"|"NcollList.dat")
        mkdir results

        # Extract the temperature and velocity profiles
        cp $musicRepoDirectory/parameters.dat .
        $musicRepoDirectory/hydro_analysis.e

        # Bin the data in NCollList.dat
        python $musicRepoDirectory/NCollHisto.py

        
        cd ../
        printf -v newDirNum "%03d" $i
        mv $dir run_$newDirNum
        i=$((i+1))
    else
        # Some of the MUSIC runs do not contain this data, delete them
        rm -r $dir
        echo "no data file, removed directory"
    fi
done