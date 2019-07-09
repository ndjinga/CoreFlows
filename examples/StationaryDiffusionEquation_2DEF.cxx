#include "StationaryDiffusionEquation.hxx"

using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 2;

	/* Mesh data */
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	int nx=5;
	int ny=5;

    /* Mesh construction */
	Mesh M(xinf,xsup,nx,yinf,ysup,ny,0); //Regular triangular mesh

	/* set the limit field for each boundary */
	double eps=1e-6;
	M.setGroupAtPlan(xsup,0,eps,"Bord1");
	M.setGroupAtPlan(xinf,0,eps,"Bord2");
	M.setGroupAtPlan(ysup,1,eps,"Bord3");
	M.setGroupAtPlan(yinf,1,eps,"Bord4");

    /* set the boundary values for each boundary */
	double T1=0;
	double T2=0;
	double T3=0;
	double T4=0;

	cout<< "Building of a regular triangular 2D mesh from a square mesh with "<< nx<<"x" <<ny<< " cells"<<endl;

    /* Create the problem */
	StationaryDiffusionEquation myProblem(spaceDim);
	myProblem.setMesh(M);

    /* set the boundary conditions */
	myProblem.setDirichletBoundaryCondition("Bord1",T1);
	myProblem.setDirichletBoundaryCondition("Bord2",T2);
	myProblem.setDirichletBoundaryCondition("Bord3",T3);
	myProblem.setDirichletBoundaryCondition("Bord4",T4);

	myProblem.setLinearSolver(GMRES,ILU);

    /* name the result file */
	string fileName = "StationnaryDiffusion_2DEF";
	myProblem.setFileName(fileName);

	/* Run the computation */
	myProblem.initialize();
	bool ok = myProblem.solveStationaryProblem();
	if (ok)
		cout << "Simulation of "<<fileName<<" is successful !" << endl;
	else
		cout << "Simulation of "<<fileName<<"  failed ! " << endl;

	cout << "------------ End of calculation !!! -----------" << endl;
	myProblem.terminate();

	return EXIT_SUCCESS;
}
