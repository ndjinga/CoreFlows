#!/usr/bin/env python
# -*-coding:utf-8 -*


import CoreFlows as cf

def SinglePhase_3DVortexTube_WithoutCone():
	spaceDim = 3;

	print( "Loading mesh of vortex tube without cone" );
	inputfile="../resources/VortexTubeWithoutCone.med";	
 
    # set the limit field for each boundary
	outletPressure =  1e5;
	inletPressure  = 10e5;
	inletTemperature  = 10e5;
	
    # physical constants
	viscosite=[0.025];

	myProblem = cf.SinglePhase(cf.Gas,cf.around1bar300K,spaceDim);
	nVar = myProblem.getNumberOfVariables();

	#Initial field creation
	print("Building initial data " ); 
	
    # Prepare for the initial condition
    
	VV_Constant = [0] * nVar

	# constant vector
	VV_Constant[0] = 1e5;
	VV_Constant[1] = 0 ;
	VV_Constant[2] = 0;
	VV_Constant[3] = 300;

    #Initial field creation
	print("Setting mesh and initial data" );
	myProblem.setInitialFieldConstant(inputfile,VV_Constant);

    # Set the boundary conditions
	myProblem.setInletPressureBoundaryCondition("Inlet flow", inletPressure, inletTemperature)
	myProblem.setOutletBoundaryCondition("Hot outlet", outletPressure)
	myProblem.setOutletBoundaryCondition("cold outlet", outletPressure)
	myProblem.setWallBoundaryCondition("Wall", inletTemperature)

    # set physical parameters
	#myProblem.setViscosity(viscosite);

    # set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
	#myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
   
    # name file save
	fileName = "3DVortexTubeWithoutCone";

    # simulation parameters
	MaxNbOfTimeStep = 3 ;
	freqSave = 1;
	cfl = 10;
	maxTime = 50;
	precision = 1e-9;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	#myProblem.setNewtonSolver(precision,20);
	#yProblem.saveConservativeField(True);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass

    # evolution
	myProblem.initialize();

	ok = myProblem.run();
	if (ok):
		print( "Simulation python " + fileName + " is successful !" );
		pass
	else:
		print( "Simulation python " + fileName + "  failed ! " );
		pass

	print( "------------ End of calculation !!! -----------" );

	myProblem.terminate();
	return ok

if __name__ == """__main__""":
    SinglePhase_3DVortexTube_WithoutCone()
