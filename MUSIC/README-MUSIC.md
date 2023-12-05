# MUSIC for use with JEWEL


In order to illustrate the use of the medium-2D code, I have collected in this directory an example MUSIC profile, along with a script that extracts the contour information and write it to files that are in the correct format for the medium-2D plugin to use.

These profiles were kindly provided by Chun Shen and are available at the following link:
https://drive.google.com/drive/folders/1rEF7Cfe2DMHmZmlUyGvMW4RjM5BijqMt

Also at that link you will find a separate README that describes what is available there and how to cite it.

## Extracting profiles from MUSIC

As you will see from the README at the link above, the hydro profiles can be extracted using a script, also provided by Chun Shen, at this link:
https://github.com/chunshen1987/hydro_analysis

Along with an additional bash script that I have written (`extractHydro.sh`), you will find the above extraction directory (`hydro_analysis`) reproduced here, but with some minor edits:
1. I have commented out lines in the `main.cpp` file that output data that is not needed by medium-2D.
2. I have edited `FluidcellStatistic.cpp` (and associated `FluidcellStatistic.h`), specifically by adding the function `outputTandUasTauvsXvsY`, modelled on the function `outputTempasTauvsX`, so that it now it outputs the temperature and velocity profiles in precisely the way read in by medium-2D.

As such, in order to run it, run the following commands
```
cd hydro_analysis
./compile.sh
./hydro_analysis.e
```

Suppose you have downloaded, from the Google Drive link above, a bunch of the runs in one of the `HYDRO_RESULTS` directories so that you have a directory that contains several folders that look like `hydro_results_##`. 
In order to extract the relevant information from these MUSIC profiles, simply copy `extractHydro.sh` into the parent directory that contains all they `hydro_results_##` folders, change the variable `musicRepoDirectory` in `extractHydro.sh` to point to the location of the `hydro_analysis` directory, and then run `extractHydro.sh`.

This will result in numerically ordered directories, containing a `results` directory with the necessary temperature and velocity profiles.
The numerical ordering is for the purposes of easily running JEWEL in parallel with numbered runs.

You can then run `jewel-2.4.0-2D` with the parameter `HYDRODIR` pointing to any of these `results` directories.

***
Note that `extractHydro.sh` is written for zsh.  If you are running bash, you might have to make minor syntactical changes.
***