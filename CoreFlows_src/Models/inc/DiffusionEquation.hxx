//============================================================================
/**
 * \file DiffusionEquation.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date 24 March 2015
 * \brief Heat diffusion equation
 * */
//============================================================================

/*! \class DiffusionEquation DiffusionEquation.hxx "DiffusionEquation.hxx"
 *  \brief Scalar heat equation for the Uranium rods temperature
 *  \details see \ref TransportEqPage for more details
 */
#ifndef DiffusionEquation_HXX_
#define DiffusionEquation_HXX_

#include "ProblemCoreFlows.hxx"

using namespace std;

class DiffusionEquation: public ProblemCoreFlows
{

public :
	/** \fn DiffusionEquation
			 * \brief Constructor for the temperature diffusion in a solid
			 * \param [in] int : space dimension
			 * \param [in] double : solid density
			 * \param [in] double : solid specific heat at constant pressure
			 * \param [in] double : solid conductivity
			 *  */

	DiffusionEquation( int dim,double rho=10000,double cp=300,double lambda=5);

	//Gestion du calcul
	void initialize();
	void terminate();//vide la mémoire et enregistre le résultat final
	bool initTimeStep(double dt);
	double computeTimeStep(bool & stop);//propose un pas de temps pour le calcul. Celà nécessite de discrétiser les opérateur (convection, diffusion, sources) et pour chacun d'employer la condition cfl. En cas de problème durant ce calcul (exemple t=tmax), renvoie stop=true
	void abortTimeStep();//efface les inconnues calculées par solveTimeStep() et reinitialise dt à 0
	bool iterateTimeStep(bool &ok);
	void save();
	void validateTimeStep();

	/** \fn setDirichletBoundaryCondition
			 * \brief adds a new boundary condition of type Dirichlet
			 * \details
			 * \param [in] string : the name of the boundary
			 * \param [in] double : the value of the temperature at the boundary
			 * \param [out] void
			 *  */
	void setDirichletBoundaryCondition(string groupName,double Temperature){
		_limitField[groupName]=LimitField(Dirichlet,-1,vector<double>(_Ndim,0),vector<double>(_Ndim,0),vector<double>(_Ndim,0),Temperature,-1,-1,-1);
	};

	void setRodDensity(double rho){
		_rho=rho;
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
protected :
	double computeDiffusionMatrix();
	double computeRHS();

	Field _fluidTemperatureField;
	bool _fluidTemperatureFieldSet, _diffusionMatrixSet;
	double _conductivity,_diffusivity, _fluidTemperature;
	double _rho;
	double _cp;
	Vector _normale;
	Matrix _DiffusionTensor;
	Vec _Tn, _deltaT, _Tk, _Tkm1, _b0;
	double _dt_diffusion, _dt_src;
};

#endif /* DiffusionEquation_HXX_ */
