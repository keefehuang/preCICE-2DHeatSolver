# preCICE-2DHeatSolver

USAGE OF 2D HEAT EQUATION SOLVER:


C++ Code:
- Contained in C++ Folder.
- Access the precice-configuration file present in \ConfigurationFiles:
	- Contains implicit and explicit configuration file for implicit and explicit coupling
- Build code using the SContruct file available
- Run using:
	 -./LeftSolver
	 -./RightSolver
- Configuration folder present in TempLeft folder as "configFile.cpp"

Python Code:
- Contained in the Python Folder
- precice-configuration file is present in the same folder:
	- Contains implicit and explicit configuration file for implicit and explicit coupling
- No building required. Activate the precice_tube environment in conda (see preCICE wiki for installation and use)
- Run using:
	- python LeftSolver.py precice-coniguration-XXXXX.xml
	- python RightSolver.py precice-coniguration-XXXXX.xml
- Configuration folder present as "configuration.py"


Configuration Variables
	- ConfigFileName: 	Name and location of configuration file
	- n_elem: 			value of N. Note that the grid is discretized equidistantly in the x- and y- axis with 2N+3 nodes 						in the x-direction and N+2 nodes in the y-direction
	- init_temp:		Initial temperature of the grid elements not at the boundary conditions
	- bound_temp:		Temperature at the dirichlet boundary conditions
	- steps:			Total number of timesteps. Note that this value is only used to keep track of the time step. The 						actual time step is determined in the preCICE condiguration file
	- dt:				Timestep used in simulation. Again note that this is only used for tracking purposes. The actual 						dt is determined in the preCICE configuration file.
	- scaling:			Arbitrary scaling factor applied to the Lagrange equation - affects BOTH solvers. Further updates 						required to manually manipulate the factor for each partitioned half (Not yet implemented for 						  C++)
	- style:			Used to switch from implict and explict time-stepping (Not yet implemented for C++)