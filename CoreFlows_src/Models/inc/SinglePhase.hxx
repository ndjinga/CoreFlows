//============================================================================
/**
 * \file DriftModel.hxx
 * \author Michael NDJINGA, Kieu Nguyen
 * \version 1.0
 * \date 01 janv. 2016
 * \brief The compressible Navier-Stokes equations
 * */
//============================================================================

/*! \class SinglePhase SinglePhase.hxx "SinglePhase.hxx"
 *  \brief The compressible Navier-Stokes equations
 *  \details The model consists in one mass, one momentum and one energy equation, see \ref NSModelsPage for more details
 */
#ifndef SINGLEPHASE_HXX_
#define SINGLEPHASE_HXX_

#include "ProblemFluid.hxx"

class SinglePhase : public ProblemFluid{
public :
	/** \fn SinglePhase
	 * \brief Constructor for the Navier-Stokes system
	 * \param [in] phaseType : \ref Liquid or \ref Gas
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 *  */
	SinglePhase(phaseType fluid, pressureEstimate pEstimate,int dim);
	//! system initialisation
	void initialize();

	//fonctions d'echange de flux
	//	void getOutputField(const Vec &Flux, const string Champ, const int numBord)=0;//, PetscInt *indices_Flux, PetscInt *indices_Bord, const long range)=0;
	//	double trace(const int &numBord, Vec &out)=0;
	void testConservation();

	void save();

	/** \fn setIntletBoundaryCondition
	 * \brief adds a new boundary condition of type Inlet
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the value of the x component of the velocity at the boundary
	 * \param [in] double : the value of the y component of the velocity at the boundary
	 * \param [in] double : the value of the z component of the velocity at the boundary
	 * \param [out] void
	 *  */
	void setInletBoundaryCondition(string groupName,double Temperature,double v_x=0, double v_y=0, double v_z=0){
		_limitField[groupName]=LimitField(Inlet,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),Temperature,-1,-1,-1);
	};
	/** \fn setIntletPressureBoundaryCondition
	 * \brief adds a new boundary condition of type InletPressure
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the pressure at the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [out] void
	 *  */
	void setInletPressureBoundaryCondition(string groupName, double pressure,double Temperature){
		_limitField[groupName]=LimitField(InletPressure,pressure,vector<double>(0,0),vector<double>(0,0),vector<double>(0,0),Temperature,-1,-1,-1);
	};
	/** \fn setWallBoundaryCondition
	 * \brief adds a new boundary condition of type Wall
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the value of the x component of the velocity at the boundary
	 * \param [in] double : the value of the y component of the velocity at the boundary
	 * \param [in] double : the value of the z component of the velocity at the boundary
	 * \param [out] void
	 *  */
	void setWallBoundaryCondition(string groupName,double Temperature,double v_x, double v_y=0, double v_z=0){
		_limitField[groupName]=LimitField(Wall,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),Temperature,-1,-1,-1);
	};

protected :
	Field _Vitesse;

	//!calcule l'etat de Roe de deux etats
	void convectionState( const long &i, const long &j, const bool &IsBord);
	//!calcule la matrice de convection de l'etat interfacial entre deux cellules voisinnes
	void convectionMatrices();
	//!Calcule le flux pour un état et une porosité et une normale donnés
	Vector convectionFlux(Vector U,Vector V, Vector normale, double porosity);
	//!calcule la matrice de diffusion de l'etat interface pour la diffusion
	void diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord);
	//!Computes the source vector associated to the cell i
	void sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i);
	//!Computes the pressure loss associated to the face ij
	void pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj);
	//!Computes the contribution of the porosity gradient associated to the face ij to the source term
	void porosityGradientSourceVector();
	//!Calcule la jacobienne de la CL convection
	void jacobian(const int &j, string nameOfGroup,double * normale);
	//!Calcule la jacobienne de la CL de diffusion
	void jacobianDiff(const int &j, string nameOfGroup);
	//!Calcule l'etat fictif a la frontiere
	void setBoundaryState(string nameOfGroup, const int &j,double *normale);// delete &nf Kieu
	//!Adds the contribution of diffusion to the RHS
	void addDiffusionToSecondMember(const int &i,const int &j,bool isBord);
	//!Computes the interfacial flux for the VFFC formulation of the staggered upwinding
	Vector staggeredVFFCFlux();
	//!Computes the matrices A^+ and A^- for the VFFC formulation of the staggered upwinding
	void staggeredVFFCMatrices(double u_n);
	//!Compute the corrected interfacial state for lowMach, pressureCorrection and staggered versions of the VFRoe formulation
	void applyVFRoeLowMachCorrections();
	//!Special preconditioner based on a matrix scaling strategy
	void computeScaling(double offset);
	//!Calcule les saut de valeurs propres pour la correction entropique
	void entropicShift(double* n);
	// Fonctions utilisant la loi d'etat
	void consToPrim(const double *Ucons, double* Vprim,double porosity=1);
	void Prim2Cons(const double *V, const int &i, double *U, const int &j);

};
#endif /* SINGLEPHASE_HXX_*/
