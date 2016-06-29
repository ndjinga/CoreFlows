#!/usr/bin/env python
# -*-coding:utf-8 -*

import CoreFlows as cf



def SinglePhase_2DSphericalExplosion_unstructured():

	print( "Loading unstructured mesh " );
	inputfile="../../CoreFlows_src/examples/resources/BoxWithMeshWithTriangularCells";
	fieldName="Initial variables for spherical explosion";
	spaceDim=2
	
	myProblem = cf.SinglePhase(cf.Liquid,cf.around155bars600K,spaceDim);
	nVar=myProblem.getNumberOfVariables();

	# Initial field creation
	print ("Setting mesh and initial data " ) ;
	myProblem.setInitialField(inputfile,fieldName,0);

	# set the boundary conditions
	wallVelocityX=0;
	wallVelocityY=0;
	wallTemperature=563;

	myProblem.setWallBoundaryCondition("GAUCHE", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("DROITE", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("HAUT", wallTemperature, wallVelocityX, wallVelocityY);
	myProblem.setWallBoundaryCondition("BAS", wallTemperature, wallVelocityX, wallVelocityY);

	# set the numerical method
	myProblem.setNumericalScheme(cf.upwind, cf.Explicit);
 	myProblem.setLinearSolver(cf.GMRES,cf.ILU,True);
	myProblem.setEntropicCorrection(False);
	myProblem.setWellBalancedCorrection(False);
    
	# name file save
	fileName = "2DSphericalExplosion_unstructred";

	# parameters calculation
	MaxNbOfTimeStep = 3 ;
	freqSave = 5;
	cfl = 0.5;
	maxTime = 5;
	precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,20);
	myProblem.saveConservativeField(True);
	if(spaceDim>1):
		myProblem.saveVelocity();
		pass


	# evolution
	myProblem.initialize();
	print("Running python "+ fileName );

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
    SinglePhase_2DSphericalExplosion_unstructured()
