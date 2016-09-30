#include "SinglePhase.hxx"
#include "DriftModel.hxx"
#include "FiveEqsTwoFluid.hxx"
#include "IsothermalTwoFluid.hxx"
#include "TransportEquation.hxx"
#include "DiffusionEquation.hxx"

#include <iostream>

using namespace std;

int main()
{
	//Preprocessing: mesh and group creation
	cout << "Building cartesian mesh" << endl;
	double xinf=0.0;
	double xsup=1.0;
	int nx=2;
	Mesh M(xinf,xsup,nx);
	double eps=1.E-8;
	M.setGroupAtPlan(xsup,0,eps,"Neumann");
	M.setGroupAtPlan(xinf,0,eps,"Neumann");
	int spaceDim = M.getSpaceDimension();

	// set the limit field for each boundary
	LimitField limitNeumann;
	limitNeumann.bcType=Neumann;
	map<string, LimitField> boundaryFields;
	boundaryFields["Neumann"] = limitNeumann;

	DriftModel  myProblem(around155bars600K,spaceDim);
	int nbPhase = myProblem.getNumberOfPhases();
	int nVar = myProblem.getNumberOfVariables();
	Field VV("Primitive", CELLS, M, nVar);//3+spaceDim*nbPhase

	// Prepare for the initial condition
	Vector VV_Left(nVar),VV_Right(nVar);
	double discontinuity = (xinf+xsup)/2.;
	VV_Left(0) = 0.5; VV_Right(0) = 0.2;
	VV_Left(1) = 155e5; VV_Right(1) = 155e5;
	for (int idim=0; idim<spaceDim;idim++){
		VV_Left(2+idim) = 1;VV_Right(2+idim) = 1;
	}
	VV_Left(2+spaceDim) = 573;
	VV_Right(2+spaceDim) = 618;

	//Initial field creation
	cout << "Building initial data" << endl;
	myProblem.setInitialFieldStepFunction(M,VV_Left,VV_Right,discontinuity);

	//set the boundary fields
	myProblem.setBoundaryFields(boundaryFields);

	// set the numerical method
	myProblem.setNumericalScheme(staggered, Implicit);
	myProblem.setWellBalancedCorrection(true);
	myProblem.setNonLinearFormulation(VFFC);

	// name file save
	string fileName = "RiemannProblem";

	//numerical parameters
	unsigned MaxNbOfTimeStep =1 ;
	int freqSave = 1;
	double cfl = 0.95;
	double maxTime = 1;
	double precision = 1e-6;

	myProblem.setCFL(cfl);
	myProblem.setPrecision(precision);
	myProblem.setMaxNbOfTimeStep(MaxNbOfTimeStep);
	myProblem.setTimeMax(maxTime);
	myProblem.setFreqSave(freqSave);
	myProblem.setFileName(fileName);
	myProblem.setNewtonSolver(precision,1);

	// set display option to monitor the calculation
	//myProblem.setVerbose( true);
	myProblem.saveConservativeField(true);

	// evolution
	myProblem.initialize();

	bool ok = myProblem.run();
	if (ok)
		cout << "Simulation "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation -----------" << endl;
	myProblem.terminate();

	return ok;
}
