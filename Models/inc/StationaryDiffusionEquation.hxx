//============================================================================
/**
 * \file StationaryDiffusionEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date June 2019
 * \brief Stationary heat diffusion equation with finite element. 
 * -\lambda\Delta T=\Phi + \lambda_{sf} (T_{fluid}-T)
 * */
//============================================================================

/*! \class StationaryDiffusionEquation StationaryDiffusionEquation.hxx "StationaryDiffusionEquation.hxx"
 *  \brief Scalar stationary heat equation for the Uranium rods temperature
 *  \details see \ref TransportEqPage for more details
 * -\lambda\Delta T=\Phi(T) + \lambda_{sf} (T_{fluid}-T)
 */
#ifndef StationaryDiffusionEquation_HXX_
#define StationaryDiffusionEquation_HXX_

#include "ProblemCoreFlows.hxx"
#include "Node.hxx"

using namespace std;

class StationaryDiffusionEquation
{

public :
	/** \fn StationaryDiffusionEquation
			 * \brief Constructor for the temperature diffusion in a solid
			 * \param [in] int : space dimension
			 * \param [in] double : solid density
			 * \param [in] double : solid specific heat at constant pressure
			 * \param [in] double : solid conductivity
			 *  */

	StationaryDiffusionEquation( int dim,bool FECalculation=true,double lambda=1);

    void setMesh(const Mesh &M);
    void setLinearSolver(linearSolver kspType, preconditioner pcType);
    void setFileName(string fileName){
	_fileName = fileName;
    }
    bool solveStationaryProblem();
    
	//Gestion du calcul
	void initialize();
	void terminate();//vide la mémoire et enregistre le résultat final
	double computeDiffusionMatrix(bool & stop);
    double computeTimeStep(bool & stop);//For coupling calculations
	bool iterateNewtonStep(bool &ok);
	void save();

    /* Boundary conditions */
	void setBoundaryFields(map<string, LimitField> boundaryFields){
		_limitField = boundaryFields;
    };
	/** \fn setDirichletBoundaryCondition
			 * \brief adds a new boundary condition of type Dirichlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setDirichletBoundaryCondition(string groupName,double Temperature){
		_limitField[groupName]=LimitField(Dirichlet,-1, vector<double>(_Ndim,0),vector<double>(_Ndim,0),
                                                        vector<double>(_Ndim,0),Temperature,-1,-1,-1);
	};

	void setConductivity(double conductivite){
		_conductivity=conductivite;
	};
	void setFluidTemperatureField(Field coupledTemperatureField){
		_fluidTemperatureField=coupledTemperatureField;
		_fluidTemperatureFieldSet=true;
	};
	void setFluidTemperature(double fluidTemperature){
	_fluidTemperature=fluidTemperature;
	}
	Field& getRodTemperatureField(){
		return _VV;
	}
	Field& getFluidTemperatureField(){
		return _fluidTemperatureField;
	}
	void setDiffusiontensor(Matrix DiffusionTensor){
		_DiffusionTensor=DiffusionTensor;
	};

	/** \fn setHeatPowerField
	 * \brief set the heat power field (variable in space)
	 * \details
	 * \param [in] Field
	 * \param [out] void
	 *  */
	void setHeatPowerField(Field heatPower){
		_heatPowerField=heatPower;
		_heatPowerFieldSet=true;
	}

	/** \fn setHeatPowerField
	 * \brief set the heat power field (variable in space)
	 * \details
	 * \param [in] string fileName (including file path)
	 * \param [in] string fieldName
	 * \param [out] void
	 *  */
	void setHeatPowerField(string fileName, string fieldName){
		_heatPowerField=Field(fileName, CELLS,fieldName);
		_heatPowerFieldSet=true;
	}

protected :
	//Main unknown field
	Field _VV;

	int _Ndim;//space dimension
	int _nVar;//Number of equations to sole
	int _Nmailles;//number of cells for FV calculation
	int _neibMaxNbCells;//maximum number of cells around a cell
	int _Nnodes;//number of nodes for FE calculation
	int _neibMaxNbNodes;//maximum number of nodes around a node
    //Data for FE calculation
	int _NinteriorNodes;//number of interior nodes for FE calculation
	int _NboundaryNodes;//number of boundary nodes for FE calculation
    bool _FECalculation;
    std::vector< int > _boundaryNodeIds;/* List of boundary nodes*/
    std::vector< int > _interiorNodeIds;/* List of interior nodes*/

	Mesh _mesh;
    bool _meshSet;
	bool _initializedMemory;
    
	double _precision;
	double _precision_Newton;
	double _erreur_rel;//norme(Uk+1-Uk)
	string _fileName;//name of the calculation
	ofstream * _runLogFile;//for creation of a log file to save the history of the simulation

	//Linear solver and petsc
	KSP _ksp;
	KSPType _ksptype;
	PC _pc;
	PCType _pctype;
	string _pc_hypre;
	int _maxPetscIts;//nombre maximum d'iteration gmres autorisé au cours d'une resolution de système lineaire
	int _maxNewtonIts;//nombre maximum d'iteration de Newton autorise au cours de la resolution d'un pas de temps
	int _PetscIts;//the number of iterations of the linear solver
	int _NEWTON_its;
	Mat  _A;//Linear system matrix
	Vec _b;//Linear system right hand side
	double _MaxIterLinearSolver;//nombre maximum d'iteration gmres obtenu au cours par les resolution de systemes lineaires au cours d'un pas de tmeps
	bool _conditionNumber;//computes an estimate of the condition number

	map<string, LimitField> _limitField;

	bool _diffusionMatrixSet;
	Vector _normale;
	Matrix _DiffusionTensor;
	Vec _deltaT, _Tk, _Tkm1, _b0;
	double _dt_src;
    
	//Heat transfert variables
	double _conductivity, _fluidTemperature;
	Field _heatPowerField, _fluidTemperatureField;
	bool _heatPowerFieldSet, _fluidTemperatureFieldSet;
	double _heatTransfertCoeff, _heatSource;

	//Display variables
	bool _verbose, _system;
    //saving parameters
	string _path;//path to execution directory used for saving results
	saveFormat _saveFormat;//file saving format : MED, VTK or CSV

	double computeRHS(bool & stop);
	double computeDiffusionMatrixFV(bool & stop);
    /* Functions for finite element method */
    Vector gradientNodal(Matrix M, vector< double > v);//gradient of nodal shape functions
	double computeDiffusionMatrixFE(bool & stop);
    int fact(int n);
};

#endif /* StationaryDiffusionEquation_HXX_ */
