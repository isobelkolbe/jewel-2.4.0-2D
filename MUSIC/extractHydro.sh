#!/bin/zsh
# file name: extractHydro.sh
# Written by Isobel Kolbe, 2023
# This script will use a slightly modified version of Chun Shen's `hydro_analysis`
# script to extract temperature and velocity profiles for use with jewel-2.4.0-2D.

# Location of the hydro analysis code that will do the extraction
musicRepoDirectory=$PWD/../hydro_analysis


setopt extended_glob 
# Loop over all the run directories
dirList=(`ls -d */`)
i=0


for dir in $dirList[@]
do 
    echo $dir
    
    dataFile=$dir/evolution_all_xyeta.dat
    if test -f "$dataFile"; then
        cd $dir 
        echo "In " $PWD 
        # Clean the directory
        rm -- *~(evolution_all_xyeta.dat|NcollList*) 
        mkdir results
        mv NcollList* NcollList.dat
        mv evolution_all_xyeta.dat results/

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