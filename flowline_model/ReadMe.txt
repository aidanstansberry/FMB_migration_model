!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Flowline Model Documentation !!!
!!! author : Aidan Stansberry    !!!
!!! August 9 2021                !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 1) How to use this model     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	a) Get ElmerFEM software https://www.csc.fi/web/elmer/binaries
	b) convert the mesh to a format usable by elmer
		$ gmsh -2 1k_mesh.geo 
		$ ElmerGrid 14 2 1k_mesh.msh
	c) Aquire necessary data sets
	d) run each step (1,2,3,4) to completion in order. .vtu output can be easily viewed in paraview, and from there data can be saved in a .csv format.
		$ ElmerSolver step.sif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!! 2) File Descriptions         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	a) 1k_mesh.geo/1k_mesh.msh
		Used to create the 2-D mesh for the flowline. A 2 layer rectangle with a top layer (ice) thickness of 20 m and a bottom layer (bedrock) thickness of 3000 m and a length of 480 km is defined. 


	b) DummySolver/DummySolver.F90
		A user function necessary for the creation and storage of several variables used in the model.

	c) smb_model_jake_data/smb_model_jake_data.F90
		user functions needed to calculate the surface mass balace for the flowline. Adapted from a surface mass balance model written in python by jake downs hence the name


	d) step1_initial_steady_state.sif
		With initial ice sheet geometry configuation run a steady state temperature iteration

	
	e) step2_run_to_steady_profile_dp_46.sif
		Allows ice sheet suface to evolve until desired initial terminus position is achieved (~3.5 kyrs)


	f) step3_initial_steady_state_temp.sif
		Get the initial internal temperature of the ice sheet and bedrock to steady state


	g) step4_run_climate_scenario.sif
		Runs the transient flowline evolution simulation from 11.4 kyr to present driven by temperature and precipitation anomalies.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!! 4) Required Data Sets        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Several Data sets were to run the model. They are described here

	a) ModernPrecip.dat
		The modern surface precipitation data along the flowline.

		first column: distance (m) | additional column(s): monthly precip (m w.e. (meters water equivalent))


	b) ModernSuface.dat
		The modern surface elevation data along the flowline.

		first column: distance (m) | additional column(s): elevation (m)


	c) ModernTemp.dat 
		Modern monthly suface temperature data along the flowline.

		first column: distance (m) | additional column(s): temperature (degC)


	d) ModernTempAvg
		Yearly average modern surface temperature data along the flowline.

		first column: distance (m) | additional column(s): temperature (degC)


	e) monthlydps.dat
		monthly precipitation anomaly data for the flowline
		
		first column: time (years) | additional column(s): monthly precip (m w.e. (meters water equivalent))

	
	f) montlydts.dat
		monthly temperature anomaly data for the flowine

		first column: time (years) | additional column(s): delta temperature (degC)

	
	g) north_initial_suface.dat
		ice sheet suface elvation to initialize the model

		first column: distance (m) | additional column(s): elevation (m)


	h) ModernPrecip.dat
		The modern surface precipitation data along the flowline.

		first column: distance (m) | additional column(s): monthly precip (m w.e. (meters water equivalent))


	i) ModernSuface.dat
		The modern surface elevation data along the flowline.

		first column: distance (m) | additional column(s): elevation (m)


	j) ModernTemp.dat 
		Modern monthly suface temperature data along the flowline.

		first column: distance (m) | additional column(s): temperature (degC)


	k) ModernTempAvg.dat
		Yearly average modern surface temperature data along the flowline.

		first column: distance (m) | additional column(s): temperature (degC)


	l) monthlydps.dat
		monthly precipitation anomaly data for the flowline
		
		first column: time (years) | additional column(s): monthly precip (m w.e. (meters water equivalent))

	
	m) montlydts.dat
		monthly temperature anomaly data for the flowine

		first column: time (years) | additional column(s): delta temperature (degC)

	
	n) north_initial_suface.dat
		ice sheet suface elvation to initialize the model

		first column: distance (m) | additional column(s): elevation (m)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!! 5) Additional Model Runs     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	a) sensitivity to bedrock topography
		Bed.dat can be changed to reflect the new bedrock configurations.

	b) sensitivity to geothermal heat flux
		This set of files can be used for these runs, but geothermal heat flux in each of the runs must be changed (e.g. 30 mW/m^2 --> 50 mW/m^2)

	c) Removal of the thermally active bedrock. Go to directory flowline_model_no_bedrock











