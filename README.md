# jewel-2.4.0-2D
A plugin for use with jewel-2.4.0 that enables the use of any (2+1)D medium profile.

## Installation
There are two simple ways to install this plugin
1. Modify an existing JEWEL installation
2. Clone this full repository

### Installation method (1) - Modify an existing JEWEL installation

The most intuitive way to include this plugin is to start with a working installation of jewel-2.4.0.  Then make the following modifications:

1. Download the file `medium-2D.f` from this repository and copy it to your jewel-2.4.0 directory (the same directory in which `medium-simple.f` is).
2. Add the following to your `Makefile` (an example Makefile is in this repository as well):
```
    jewel-2.4.0-2D: medium-2D.o jewel-2.4.0.o pythia6425mod-lhapdf6.o meix.o
        $(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++
```       
Don't forget to add `jewel-2.4.0-2D` to the `all:` line.
3. I also suggest making the following minor modification to the `jewel-2.4.0.f` code.  Depending on the naming structure of your medium profiles, it might not matter, but if you will be running on a cluster or have any other reason for needing larger file names or absolute file locations, modify all instances of `FILEMED`, `buffer`, and `value` to be of type `CHARACTER*300 (or whatever length ensures that the entire name fits into the variable).

You should now be able to do `make`, which should result in the usual `jewel-2.4.0-vac` and `jewel-2.4.0-simple` programs, as well as the new `jewel-2.4.0-2D` program.

### Installation method (2) - Clone this repository

Alternatively, if you have a working LHAPDF-6 installation, you could clone this entire
repository, edit the LHAPDF line in the Makefile as per the JEWEL installation instructions,
and build.

More concretely, follow these steps:
1. Clone this repository by doing
```
git clone https://github.com/isobelkolbe/jewel-2.4.0-2D
```
2. Find out the exact location of your LHAPDF6 installation and edit the `LHAPDF_PATH` variable in the `Makefile` to point to the LHAPDF library, i.e. the path to `LHAPDF-6.#.#/lib`.

## Running 

The 2D medium profiles are read in by the subroutines `readtemps(filename)` and `redvelocities(filename)`.  You are, of course, welcome to edit those routines to fit the format of your medium profiles, but it is probably easier to use a bash script to change the format of your profiles to match what the subroutines expect.  For example, in the repository, you will find a script that modifies publicly available MUSIC profiles for use with `jewel-2.4.0-2D`.

The subroutines expect the input data to have the following properties:
1. Temperature and velocity information seperately. They do not have to be on the same (x,y) grid.
2. Temperature (velocity) files in the format `Tcontour##.###` (`Vcontour##.###`), where `##.###` is the time stamp (in fm/c)
3. Each temperatur file should contain space or comma separated columns.  The first coloumn contains x-values, the second y-values, and the third temperature value at that (x,y)-position. For the velocity files, the third and fourth columns should contain the x-component and y-component of the fluid velocity.
4. The files may contain additional columns, they will not be read in.
5. The x- and y- values should be uniformly separated, increasing, and the same in all files.
6. The first line is assumed to be a column heading and is not read in.

If available (highly recommended for small systems), it is possible also to include a probability distribution of the binary collisions (for example those used to source hydrodynamical simmulations).  If providing such a probability distribution, `jewel-2.4.0-2D` will use it to sample the position of the initial hard scattering (instead of a Glauber overlap function).  The name of the file containing this data can be passed as a parameter in the medium parameter file, and should have the same format as the temperature and velocity profiles: space or comma separated columns containg the x and y positions with the probability density at that location (an integer value).  Again, an script that extracts this information from MUSIC data is provided.
