--- Contents ---

1: Introduction

2: Citation Details

3: Usage
   3.1: Dependencies
   3.2: Compilation
   3.3: Output Files

4: Correspondence

5: License


--- 1: Introduction ---

Critical Configurations of Hard Disks on the Torus is a collection of 
scripts and C++ code to identify the critical configurations of hard disks 
on the flat torus. The current version contains source code to handle 
anywhere from one to twelve disks, and is parallelized using OpenMP for 
shared memory machines.

The base algorithm randomly selects the initial locations of the hard disks
on the sixty degree flat torus. The interaction potential of a pair of disks
is given by

\exp(-w \cdot (R - r))

where $R$ is the pairwise distance and $r$ is the minimum pairwise distance 
over all disks. This potential approaches the hard disk potential as $w$
increases to infinity, a feature that considerably improves the numerical
performance.

The potential energy of the system is given by summing the interaction 
potential of all pairs of disks, and the fictive energy is given by the
norm of the gradient of the potential energy. The locations of the disks
evolve by a conjugate gradient descent on the fictive energy surface,
meaning that the critical points of the potential energy surface result in
stable attractors.

When the conjugate gradient descent finishes, $w$ is increased, shifting
the positions of the critical points of the potential energy surface and
requiring the descent to be repeated. Provided that $w$ is increased 
slowly, the system is expected to remain in the same fictive energy basin.
This is intended as a workaround to avoid the difficulty with following
a gradient descent on the hard disk potential directly.

This program is open-source software, and is distributed under the GPLv2. 
It was completed on 2013/04/16 at Lawrence Livermore National Laboratory by 
Jeremy Mason, reachable at jkylemason@gmail.com.


--- 2: Citation Details ---

The algorithm and implementation of this code is not currently documented
anywhere outside of these files. If this code is used for research 
purposes or as part of other software, please acknowledge the author.


--- 3: Usage ---

The makefile generates a disk_search executable that is run from the 
command line as

$ ./disk_search [OPTIONS FILE]

An options file must be provided. A user-defined options file should be 
based on options.txt, and generally contains the following keywords and 
values:

NDisks [REQUIRED]
  A positive integer in interval [1, 12] - number of disks on the torus
NTrials [REQUIRED]
  Any positive integer - number of initial configurations attempted

NThreads [REQUIRED]
  Any positive integer - number of threads when parallelization is enabled
WIncrements [REQUIRED]
  Any positive integer - number of $w$ increments used to harden potential
EThreshold [REQUIRED]
  Any positive float - energy threshold below which configuration is 
    considered to be converged

TimeInterval [REQUIRED]
  Any positive integer - number of seconds between saves to disk
BackupInterval [REQUIRED]
  Any positive integer - number of seconds between generation of backups

The number of disks is hardcoded into the source, and must be consistent
with the value of NDisks in the options file.

If the executable finds critical_point.txt and occurrences.txt files in the
current directory, and these are for a number of disks that is compatible
with the options file, these configurations are incorporated into the list
for the current simulation. Otherwise, empty critical_point.txt and 
occurrences.txt files are generated. These are updated periodically as
the simulation proceeds to reduce disk access overhead.

The user is expected to manage archiving of the critical_point.txt and 
occurrences.txt files to reduce the chance of earlier simulation results
being overwritten. Directories are provided in results/ for this purpose.


--- 3.1: Dependencies ---

The current version of the code relies on several external libraries. 

nauty is a program for computing automorphism groups of graphs available at
http://cs.anu.edu.au/~bdm/nauty/ and should be installed at src/libraries/. 
This is required to find the canonically labelled contact graphs of the 
critical configurations. Tested using the 2.5 version.

armadillo is a comprehensive C++ linear algebra library available at 
http://arma.sourceforge.net/ and should be installed per the instructions in
the README. This is required for all internal linear algebra operations. Tested
using the 7.600.2 version.

g++ is recommended for compilation of the C++ code. The parallel version
of the simulation depends on the availability of OpenMP.

This is sufficient to compile the executable and perform the simulation.
The distribution includes files that extend the functionality further,
provided that further dependencies are available.

MATLAB is required to view the resulting configurations with the included 
script plotting.m. The MATLAB version provided in matlab/ is a prototype
of the C++ version, and is useful to set the numerical parameters of the
simulation. Finally, the MATLAB Symbolic Math Toolbox is used by 
setup/generate_functions.m to generate the source code required to run the 
simulation for different numbers of disks. Tested using the 2018a version.

Python is used to run cleanup.py and convert.py in setup/ to convert the 
files generated by generate_functions.m into valid MATLAB functions and 
C++ code. This is only required to generate the source code required to 
run the simulation for more than twelve disks. Tested using the 3.6
version.


--- 3.2: Compilation ---

The user is expected to compile the executable directly on the local 
architecture. The recommended procedure is to use the setup.sh bash script

$ ./setup.sh [NUMBER OF DISKS]

where the number of disks is an integer in the interval [1, 12]. This 
copies the appropriate source files from /setup to /src, changes NDisks in
options.txt to the appropriate value, and compiles the executable 
disk_search with parallelization enabled in the current directory.

Other compilation options are provided in the makefile. Generally, these
are used only during code development.

The user does need to modify the makefile to reflect the availability of
linear algebra libraries linked to by Armadillo. This is usually as simple
as commenting out all but one of the included LIBS definitions before 
running setup.sh.

If the user desires to run the MATLAB version provided in matlab/, a mex
function must first be compiled on the local architecture. This is done by
navigating to the matlab/ directory within MATLAB and entering

>> mex -I../src/libraries/nauty/ canon_label.c ../src/libraries/nauty/nauty.c ../src/libraries/nauty/nautil.c ../src/libraries/nauty/naugraph.c ../src/libraries/nauty/schreier.c ../src/libraries/nauty/naurng.c

to generate a MATLAB interface to nauty. The MATLAB version is run by
executing the script matlab/main.m, and provides a convenient interface to
observe the effect of changing simulation parameters.


--- 3.3: Output Files ---

The simulation periodically saves the available critical point information 
in two text files.

critical_point.txt - this file lists the known critical points in order of
  their identification. The first line is the number of disks for all of
  the critical points. The following information is then given for every
  critical point:
    the adjacency matrix for the canonically labelled contact graph
    the x and y coordinates of the disks
    the radius of all of the disks
    the index of the critical point, estimated from the number of negative 
      eigenvalues of the Hessian
  This file is appended with entries for additional critical points as they
  are identified during the simulation. The preferred means of viewing the
  resulting critical configurations is by executing the plotting.m script
  in the same directory as a critical_point.txt file and following the 
  instructions at the MATLAB prompt.

occurrences.txt - this file lists the number of observations of each of the 
  critical points listed in critial_point.txt, and is intended to give a 
  sense of when the list of critical points is nearing completion. This is
  expected to occur when all of the values in occurrences.txt are above 
  some predetermined value (~ 10) assuming that there is some continuity in 
  the rates of occurrence for the entire population of critical points.


--- 4: Correspondence ---

Please direct all inquiries to Jeremy K. Mason at jkylemason@gmail.com.


--- 5: License ---

This program is open-source software, and is distributed under the GPLv2.
Please refer to the "LICENSE.txt" file for more complete license
information.

