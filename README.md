# SpringRegister ###############################################################

## Authors #####################################################################

Primary Author: [Noah C. Benson](mailto:n@nben.net)
Principle Investigator: [Geoffrey K. Aguirre](mailto:aguirreg@mail.med.upenn.edu)

## Citation ####################################################################

This repository is part of the paper: Benson NC, Butt OH, Brainard DH,
Aguirre GK (2014) Correction of distortion in flattened
representations of the cortical surface allows prediction of V1-V3
functional organization from anatomy. PLoS Comput Biol. Submitted.

## Description #################################################################

This repository contains the C++ source code used to perform spring
registration of a flattened cortical surface to an ideal 2D model in
Benson, Butt, Brainard, Aguirre (2014).

This src directory contains the source code as well as a Makefile,
which will build the springs executable file. This program is designed
to compile and run on unix-like systems, and was tested on Mac OS X
and Ubuntu. GNU Make version 3.81 and G++ (GCC) version 4.2.1 were
used in building the executable. The standard math library (-lm) and
POSIX (-lpthread) are also required for building.

This program is intended for use in conjunction with the jobs written
by the supplementary Mathematica notebook included in the publication
cited above. The springs program is, in fact, a low-level program
designed to increase speed and parallelization of the numerical
integration used to simulate and minimize the spring system described
in Benson et al.

Many command line options to the springs program exist; these can be
examined by calling springs -h or springs --help. If you are running
simulations whose jobs were exported by the Mathematica notebook, then
these options will be specified for you, and you may simply execute
the run.sh file that is written. Note that this script expects a
certain directory structure in which <root>/jobs/<job>/run.sh is
called from <root> and in which <root>/src contains the springs
executable. This repository exemplifies this organization, and the
jobs and results directories analyzed in Benson et al. are included.

## License #####################################################################

See LICENSE file.

