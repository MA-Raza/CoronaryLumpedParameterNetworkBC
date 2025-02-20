- A working installation of solids4Foam-v2.1 with foam-extend-4.1
is required to compile the BCs and run the tuturial case.

- LumpedParameterNetworkBCTestCase is the parent folder containing
everything inside it.

- A master directory (LPNBCTestCase) containing master dictionary 
templates for all the simulation cases is given within the parent
folder.

- A directory (lumpedParameterNetworkBC) containing the codes for
developed BCs is also placed within the parent folder.

- MATLABPostprocessing directory contains the MATLAB script to
post process the results of simulations.

- To run the simulations do the following:

(1) Run the "CreateDirStruc" bash script in the parent folder, which
will create the basic directory structure for all the simulation
cases (11 in total).

(2) Run the main "Allrun" bash script in the parent folder, which will
execute "Allwmake" script in the BC directory to compile the BCs, and
then will run all the cases one after another using 8 processors in
parallel (by executing local "Allrun" scripts in each child case directory).

(3) After all the simulation cases are completed, run the MATLAB script
placed inside MATLABPostprocessing directory to postprocess the results
and draw figures.

(4) All the simulation cases and results can be cleaned by execting main
"Allclean" bash script.

(5) The directory structure can be deleted by executing the the main
"DeleteDirStruc" bash script.