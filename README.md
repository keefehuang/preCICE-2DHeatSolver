# preCICE-2DHeatSolver
Developed based on preCICE library and elastictube1d example, available: https://github.com/precice, 
for seminar: Partitioned Fluid-Structure Interaction - Summer 18

LeftSolver code for both C++ and Python implementations are commented.

USAGE OF 2D HEAT EQUATION SOLVER:

C++ Code:
- Contained in C++ Folder.
- Access the precice-configuration file present in \ConfigurationFiles:
	- Contains implicit and explicit configuration file for implicit and explicit coupling
- Build code using the SContruct file available:
	- run "scons parallel=no petsc=no"
- Run using:
	 -run "./LeftSolver"
	 -run "./RightSolver"
- Configuration folder present in TempLeft folder as "configFile.cpp"

Python Code:
- Contained in the Python Folder
- precice-configuration file is present in the same folder:
	- Contains implicit and explicit configuration file for implicit and explicit coupling
- No building required. Activate the precice_tube environment in conda (see preCICE wiki for installation and use)
- Run using:
	- python LeftSolver.py precice-coniguration-[implicit/explicit].xml
	- python RightSolver.py precice-coniguration-[implicit/explicit].xml
- Configuration folder present as "configuration.py"


Configuration Variables in configuration files:
	- ConfigFileName: 	Name and location of configuration file
	- n_elem: 			value of N. Note that the grid is discretized equidistantly in the x- and y- axis with 2N+3 nodes 						in the x-direction and N+2 nodes in the y-direction
	- init_temp:		Initial temperature of the grid elements not at the boundary conditions
	- bound_temp:		Temperature at the dirichlet boundary conditions
	- steps:			Total number of timesteps. Note that this value is only used to keep track of the time step. The 						actual time step is determined in the preCICE condiguration file
	- dt:				Timestep used in simulation. Again note that this is only used for tracking purposes. The actual 						dt is determined in the preCICE configuration file.
	- scaling:			Arbitrary scaling factor applied to the Lagrange equation - affects BOTH solvers. Further updates 						required to manually manipulate the factor for each partitioned half (Not yet implemented for 						  C++)
	- style:			Used to switch from implict and explict time-stepping (Not yet implemented for C++)


References
[1] M. Mehl, B. Uekermann, H. Bijl, D. Blom, B. Gatzhammer, and A. van Zuijlen. Parallel coupling numerics for partitioned fluid-structure interaction simulations. CAMWA, 2016.