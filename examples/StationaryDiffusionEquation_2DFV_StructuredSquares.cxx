#include "StationaryDiffusionEquation.hxx"
#include "math.h"

double pi = M_PI;
using namespace std;

int main(int argc, char** argv)
{
	int spaceDim = 2;

	/* Mesh data */
	double xinf=0.0;
	double xsup=1.0;
	double yinf=0.0;
	double ysup=1.0;
	int nx=30;
	int ny=30;

    /* Mesh construction */
	Mesh M(xinf,xsup,nx,yinf,ysup,ny); //Regular square mesh

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

	cout<< "Building of a regular square 2D mesh with "<< nx<<"x" <<ny<< " cells"<<endl;

    /* Create the problem */
    bool FEComputation=false;
	StationaryDiffusionEquation myProblem(spaceDim,FEComputation);
	myProblem.setMesh(M);

    /* set the boundary conditions */
	myProblem.setDirichletBoundaryCondition("Bord1",T1);
	myProblem.setDirichletBoundaryCondition("Bord2",T2);
	myProblem.setDirichletBoundaryCondition("Bord3",T3);
	myProblem.setDirichletBoundaryCondition("Bord4",T4);

	/* Set the right hand side function*/
	Field my_RHSfield("RHS_field", CELLS, M, 1);
    Cell Ci; 
    double x, y;
	for(int i=0; i< M.getNumberOfCells(); i++)
    {
		Ci= M.getCell(i);
		x = Ci.x();
		y = Ci.y();

		my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y);//mettre la fonction definie au second membre de l'edp
	}
	myProblem.setHeatPowerField(my_RHSfield);
	myProblem.setLinearSolver(GMRES,ILU);

    /* name the result file */
	string fileName = "StationnaryDiffusion_2DFV_StructuredSquares";
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
