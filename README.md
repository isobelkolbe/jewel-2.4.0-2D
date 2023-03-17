# jewel-2.4.0-2D
A plugin for use with jewel-2.4.0 that enables the use of any (2+1)D medium profile.


There are two simple ways to install this plugin
    (1) Modify an existing JEWEL installation
    (2) Clone this full repository

METHOD (1) - Modify an existing JEWEL installation
---------------------------------------------------
The most intuitive way to include this plugin is to start with a working installation of
jewel-2.4.0.  Then make the following modifications:
(1.1)- Download the file medium-2D.f from this repository and copy it to your jewel-2.4.0
    directory (the same directory in which medium-simple.f is).
(1.2)- Add the following to your Makefile (an example Makefile is in this repository as well):

    jewel-2.4.0-2D: medium-2D.o jewel-2.4.0.o pythia6425mod-lhapdf6.o meix.o
        $(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++
       
       Don't forget to add jewel-2.4.0-2D to the ``all:'' line.

It is then simply a matter of copying the medium-2D.f file into the 
jewel-2.4.0 directory, adding a line to the Makefile, and rebuilding.  Depending on how
you will be running JEWEL, I also recommend chan

Alternatively, if you have a working LHAPDF-6 installation, you could clone this entire
repository, edit the LHAPDF line in the Makefile as per the JEWEL installation instructions,
and build.