# jewel-2.4.0-2D
A plugin for use with jewel-2.4.0 that enables the use of any (2+1)D medium profile.

## Licensing and citations

In this repository you will find the source code for JEWEL, a plugin for JEWEL called medium-2D, example publicly available MUSIC profiles, and some useful scripts.

The plugin

If you use this code, please cite as follows:
1. For the standard `jewel-2.4.0-vac` and `jewel-2.4.0-simple` code, as per the JEWEL licence (and `README-jewel` here), cite the manual for the code, Eur.Phys.J. C74 (2014) no.2, 2762 [arXiv:1311.0048] and JHEP 1303 (2013) 080 [arXiv:1212.1599], and optionally EPJC C60 (2009) 617 [arXiv:0804.3568] for the physics.
2. For the `jewel-2.4.0-2D` plugin, [arxXiv:2303.14166]
3. If you use the MUSIC profiles that are here, the inSpire handle is Schenke:2020mbo.

## Installation
There are two simple ways to install this plugin
1. Modify an existing JEWEL installation
2. Clone this full repository

They are really equivalent and should result in the exact same setup.

### Installation method (1) - Modify an existing JEWEL installation

The most intuitive way to include this plugin is to start with a working installation of jewel-2.4.0.  Then make the following modifications:

1. Download the file `medium-2D.f` from this repository and copy it to your jewel-2.4.0 directory (the same directory in which `medium-simple.f` is).
2. Add the following to your `Makefile` (an example Makefile is in this repository as well):
    ```
        jewel-2.4.0-2D: medium-2D.o jewel-2.4.0.o pythia6425mod-lhapdf6.o meix.o
            $(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++
    ```  
    Don't forget to add `jewel-2.4.0-2D` to the `all:` line.
3. I also suggest making the following minor modification to the `jewel-2.4.0.f` code:  Depending on the naming structure of your medium profiles, it might not matter, but if you will be running on a cluster or have any other reason for needing larger file names or absolute file locations, modify all instances of `FILEMED`, `buffer`, and `value` to be of type `CHARACTER*300 (or whatever length ensures that the entire name fits into the variable).

You should now be able to do `make`, which should result in the usual `jewel-2.4.0-vac` and `jewel-2.4.0-simple` programs, as well as the new `jewel-2.4.0-2D` program.

### Installation method (2) - Clone this repository

Alternatively, if you have a working LHAPDF-6 installation, you could clone this entire
repository, edit the LHAPDF line in the Makefile as per the JEWEL installation instructions,
and build.  
This repository really just contains the JEWEL source code as it arrives downloaded from https://jewel.hepforge.org/ and modified as in method (1).

More concretely, follow these steps:
1. Clone this repository by doing
```
git clone https://github.com/isobelkolbe/jewel-2.4.0-2D
```
2. Find out the exact location of your LHAPDF6 installation and edit the `LHAPDF_PATH` variable in the `Makefile` to point to the LHAPDF library, i.e. the path to `LHAPDF-6.#.#/lib`.

You should now be able to do `make`, which should result in the usual `jewel-2.4.0-vac` and `jewel-2.4.0-simple` programs, as well as the new `jewel-2.4.0-2D` program.

## Running 

The 2D medium profiles are read in by the subroutines `readtemps(filename)` and `redvelocities(filename)`.  You are, of course, welcome to edit those routines to fit the format of your medium profiles, but it is probably easier to use a bash script to change the format of your profiles to match what the subroutines expect.  For example, in the repository, you will find a script that modifies publicly available MUSIC profiles for use with `jewel-2.4.0-2D`.

The subroutines expect the input data to have the following properties:
1. Temperature and velocity information seperately, in contours for each time-stamp.
2. Temperature (velocity) files in the format `Tcontour##.###` (`Vcontour##.###`), where `##.###` is the time stamp (in fm/c). (The routine will extract a variable of type float with 3 decimal places.)
3. Each temperature file should contain space or comma separated columns.  The first coloumn contains x-values, the second y-values, and the third the temperature value at that (x,y)-position. For the velocity files, the third and fourth columns should contain the x-component and y-component of the fluid velocity.
4. The files may contain additional columns, they will not be read in.
5. The x- and y- values should be uniformly separated, increasing, and the same in all files.
6. The first line is assumed to be a column heading and is not read in.

If available (highly recommended for small systems), it is possible also to include a probability distribution of the binary collisions (for example those used to source hydrodynamical simmulations).  If providing such a probability distribution (in the form of binned counts of N_coll), `jewel-2.4.0-2D` will use it to sample the position of the initial hard scattering (instead of a Glauber overlap function).  The name of the file containing this data can be passed as a parameter (`NCOLLHISTO`) in the medium parameter file, and should have the same format as the temperature and velocity profiles: space or comma separated columns containg the x- and y- positions in the first and second columns, and N_coll in a bin with (x,y) at the center of the bin in the third column (an integer value).  Again, an example script that extracts this information from MUSIC data is provided.
*** NB!! If you provide an `NCOLLHISTO`, it is imperative that you also provide the centrality and impact parameter boundaries associated with your medium model.  If your medium model uses a Glauber initial condition, then 1710.07098 might prove a useful resource for determining centrality and impact parameter limits. ***

The temperature and velocity contours can either be in the same directory as `jewel-2.4.0-2D`, or they can be in a different directory which is passed as a parameter (`HYDRODIR`) to the medium parameter file.
In this way it is possible to mimic event-by-event runs by having several different profiles in different directories, running some number of events on each profile, and merging the resultant outputs.

When running `jewel-2.4.0-2D` it will output a test message to ensure that the reading of the contours has gone smoothly. The quoted indices refer to the arrays stored by the program, so that, for instance, the temperatrue value at position `(1,1,1)` is not at time 1 fm and (x,y) = (1,1)fm, but rather the very first entry, so it will be referring to the top-left most temperature value at the first time stamp (so it is likely to be zero).  However, the x-value should not be zero.  The program will also throw errors if it was unable to read the input files, but it is best to make sure that you get something you expect.

Beyond these considerations the program works precisely like the rest of JEWEL, which is to say that you can now run any of the three programs with a parameter file:
```
./jewel-2.4.0-2D params/params-2D.dat
```
