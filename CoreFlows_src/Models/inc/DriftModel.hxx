//============================================================================
/**
 * \file DriftModel.hxx
 * \author Michael NDJINGA
 * \version 1.0
 * \date 01 janv. 2016
 * \brief Four equation two phase flow drift model
 * */
//============================================================================

/*! \class DriftModel DriftModel.hxx "DriftModel.hxx"
 *  \brief Four equation two phase flow drift model
 *  \details One total mass equation, one vapour mass equation, one total momentum equation, one total energy equation, see \ref DriftModelPage for more details
 */
#ifndef DRIFTMODEL_HXX_
#define DRIFTMODEL_HXX_

#include "ProblemFluid.hxx"

class DriftModel : public ProblemFluid{
public :
	/** \fn DriftModel
	 * \brief Constructor for the four equation drift model
	 * \param [in] pressureEstimate : \ref around1bar or \ref around155bars
	 * \param [in] int : mesh dimension
	 * \param [in] bool : There are two possible equations of state for each phase
	 *  */
	DriftModel( pressureEstimate pEstimate, int dim, bool _useDellacherieEOS=true);
	//! system initialisation
	void initialize();

	//fonctions d'echange de flux
	//	void getOutputField(const Vec &Flux, const string Champ, const int numBord)=0;//, PetscInt *indices_Flux, PetscInt *indices_Bord, const long range)=0;
	//	double trace(const int &numBord, Vec &out)=0;

	void testConservation();
	void save();

	/** \fn saveVoidFraction
	 * \brief saves the void fraction field in a separate file
	 * @param boolean saveField
	 * */
	void saveVoidFraction(bool saveField=true){
		_saveVoidFraction=saveField;
	}

	/** \fn saveEnthalpy
	 * \brief saves the enthalpy field in a separate file
	 * @param boolean saveEnthalpy
	 * */
	void saveEnthalpy(bool saveEnthalpy=true){
		_saveEnthalpy=saveEnthalpy;
	}

	// Boundary conditions
	/** \fn setIntletBoundaryCondition
	 * \brief adds a new boundary condition of type Inlet
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the value of the concentration at the boundary
	 * \param [in] double : the value of the x component of the velocity at the boundary
	 * \param [in] double : the value of the y component of the velocity at the boundary
	 * \param [in] double : the value of the z component of the velocity at the boundary
	 * \param [out] void
	 *  */
	void setInletBoundaryCondition(string groupName,double Temperature,double concentration, double v_x=0, double v_y=0, double v_z=0){
		_limitField[groupName]=LimitField(Inlet,-1,vector<double>(1,v_x),vector<double>(1,v_y),vector<double>(1,v_z),Temperature,-1,-1,concentration);
	};
	/** \fn setIntletPressureBoundaryCondition
	 * \brief adds a new boundary condition of type InletPressure
	 * \details
	 * \param [in] string : the name of the boundary
	 * \param [in] double : the value of the pressure at the boundary
	 * \param [in] double : the value of the temperature at the boundary
	 * \param [in] double : the value of the concentration at the boundary
	 * \param [out] void
	 *  */
	void setInletPressureBoundaryCondition(string groupName, double pressure,double Temperature,double concentration){
		_limitField[groupName]=LimitField(InletPressure,pressure,vector<double>(0,0),vector<double>(0,0),vector<double>(0,0),Temperature,-1,-1,concentration);
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

	/** \fn computeNewtonVariation
	 * \brief Builds and solves the linear system to obtain the variation Vkp1-Vk in a Newton scheme using primitive variables
	 * @param
	 * */
	void computeNewtonVariation();

	/** \fn iterateTimeStep
	 * \brief calls computeNewtonVariation to perform one Newton iteration and tests the convergence
	 * @param
	 * @return boolean ok is true is the newton iteration gave a physically acceptable result
	 * */
	bool iterateTimeStep(bool &ok);

protected :
	double _khi, _ksi, _kappa;//mixture pressure derivatives with regard to rhom, m_c and rhom_em
	Field _Vitesse, _VoidFraction, _Enthalpy;
	bool _saveVoidFraction,_saveEnthalpy;
	/** \fn convectionState
	 * \brief calcule l'etat de Roe de deux etats
	 * @param i,j sont des entiers qui correspondent aux numeros des cellules à gauche et à droite de l'interface
	 * @param IsBord est un booléen qui dit si la cellule i est sur le bord
	 * @return  l'état de Roe de i et j*/
	void convectionState( const long &i, const long &j, const bool &IsBord);

	/** \fn convectionMatrices
	 * \brief calcule la matrice de convection de l'etat interfacial entre deux cellules voisinnes */
	void convectionMatrices();

	/** \fn Flux
	 * \brief Computes the convection flux F(U) projected on a vector n
	 * @param U is the conservative variable vector (total density, vapour partial density, total momentum, total energy)
	 * @param V is the extended primitive variable vector (vapour mass concentration, pressure, mixture velocity, temperature, vapour volume fraction, relative velocity, vapour density, liquid density mixture enthalpy)
	 * @param normal is a unit vector giving the direction where the convection flux matrix F(U) is projected
	 * @param porosity is the ration of the volume occupied by the fluid in the cell (default value is 1)
	 * @return The convection flux projected in the direction given by the normal vector: F(U)*normal */
	Vector convectionFlux(Vector U,Vector V, Vector normale, double porosity);

	/** \fn diffusionStateAndMatrices
	 * \brief calcule la matrice de diffusion de l'etat interface pour la diffusion
	 * @param i,j sont des entiers qui correspondent aux indices  des cellules de gauche et droite respectivement
	 * @param IsBord: bollean telling if (i,j) is a boundary face
	 * @return la matrice de diffusion de l'etat interface pour la diffusion calculé */
	void diffusionStateAndMatrices(const long &i,const long &j, const bool &IsBord);

	/** \fn sourceVector
	 * \brief Ajoute au second membre la contribution de la gravité
	 * @param Si vecteur de double correspond au Snd membre
	 * @param Ui vecteur de double contient les variables conservatives
	 * @param Vi vecteur de double contient les variables primitives
	 * @param i int representing the cell number. Used for reading power field
	 * @return Snd membre mis à jour en ajoutant la contribution de la gravité */
	void sourceVector(PetscScalar * Si,PetscScalar * Ui,PetscScalar * Vi, int i);

	//Computes the pressure loss associated to the face ij
	/** \fn pressureLossVector
	 * \brief Computes the contribution of pressure loss terms in the source term
	 * \Details pure virtual function, overloaded by each model
	 * @param pressureLoss output vector containing the pressure loss contributions
	 * @param K, input pressure loss coefficient
	 * @param Ui input primitive vectors
	 * @param Vi input conservative vectors
	 * @param Uj input primitive vectors
	 * @param Vj input conservative vectors
	 * @return
	 * */
	void pressureLossVector(PetscScalar * pressureLoss, double K, PetscScalar * Ui, PetscScalar * Vi, PetscScalar * Uj, PetscScalar * Vj);

	//Computes the contribution of the porosity gradient associated to the face ij to the source term
	/** \fn porosityGradientSourceVector
	 * \brief Computes the contribution of the porosity variation in the source term
	 * \Details pure virtual function, overloaded by each model
	 * @param porosityGradientVector output vector containing the porosity variation contributions to the source term
	 * @param Ui input primitive vectors celli
	 * @param Vi input conservative vectors celli
	 * @param porosityi input porosity value celli
	 * @param deltaxi input diameter celli
	 * @param Uj input primitive vectors cellj
	 * @param Vj input conservative vectors cellj
	 * @param porosityj input porosity value cellj
	 * @param deltaxj input diameter cellj
	 * @return
	 * */
	void porosityGradientSourceVector();

	/* \fn gravityMatrix
	 * \brief matrice de gravite
	 * @param
	 * @return
	 void gravityMatrix(); */

	/** \fn jacobian
	 * \brief Calcule la jacobienne de la CL convection
	 * @param j est un entier constant correpond  à l'indice de la cellule sur le bord
	 * @param nameOfGroup est une chaine de caractère qui correspond au nom du bord
	 * 		 * @return Calcule la jacobienne de la CL convection calculé*/
	void jacobian(const int &j, string nameOfGroup,double * normale);

	/** \fn jacobianDiff
	 * \brief Calcule la jacobienne de la CL de diffusion
	 * @param j est un entier constant correpond  à l'indice de la cellule sur le bord
	 * @param nameOfGroup est une chaine de caractère qui correspond au nom du bord
	 * @return la jacobienne de la CL de diffusion calculé */
	void jacobianDiff(const int &j, string nameOfGroup);

	/** \fn setBoundaryState
	 * \brief Calcule l'etat fictif a la frontière
	 * @param j est un entier constant correpond  à l'indice de la cellule sur le bord
	 * @param nameOfGroup est une chaine de caractère qui correspond au nom du bord
	 * @return l'etat fictif a la frontière calculé */
	void setBoundaryState(string nameOfGroup, const int &j,double * normale);// delete &nf Kieu

	/** \fn addDiffusionToSecondMember
	 * \brief Compute the contribution of the diffusion operator to the right hand side of the system
	 * \Details this function is pure virtual, and overloaded in each physical model class
	 * @param i left cell number
	 * @param j right cell number
	 * @param boolean isBoundary is true for a boundary face (i,j) and false otherwise
	 * */
	void addDiffusionToSecondMember(const int &i,const int &j,bool isBoundary);
	//!Computes the interfacial flux for the VFFC formulation of the staggered upwinding
	Vector staggeredVFFCFlux();
	//!Computes the matrices A^+ and A^- for the VFFC formulation of the staggered upwinding
	void staggeredVFFCMatrices(double u_mn);
	//!Computes the matrices A^+ and A^- for the VFFC formulation of the staggered upwinding using Primitive Variables
	void staggeredVFFCMatricesPrimitiveVariables(double u_mn);
	//!Compute the corrected interfacial state for lowMach, pressureCorrection and staggered versions of the VFRoe formulation
	void applyVFRoeLowMachCorrections();
	//!Calcule les saut de valeurs propres pour la correction entropique
	void entropicShift(double* n);

	/** \fn computeScaling
	 * \brief Special preconditioner based on a matrix scaling strategy
	 * @param offset est un double, correspond au la plus grande valeure propre "Maxvp"
	 */
	void computeScaling(double offset);

	// Fonctions utilisant la loi d'etat 

	/** \fn consToPrim
	 * \brief computes the primitive vector state from a conservative vector state
	 * @param Ucons : conservative variable vector
	 * @pram Vprim : primitive variable vector
	 * @param porosity is the porosity coefficient in case of a porous modeling
	 * */
	void consToPrim(const double *U, double* V,double porosity=1);

	/** \fn primToCons
	 * \brief computes the conservative vector state from a primitive vector state
	 * @param U : conservative variable vector, may contain several states
	 * @pram V : primitive variable vector, may contain several states
	 * @param i : index of the conservative state in the vector U
	 * @param j : index of the primitive state in the vector V
	 * 	 */
	void primToCons(const double *V, const int &i, double *U, const int &j);

	/** \fn primToConsJacobianMatrix
	 * \brief computes the jacobian matrix of the cons->prim function
	 * @pram V : primitive vector state
	 * 	 */
	void primToConsJacobianMatrix(double *V);

	/** \fn computeExtendedPrimState
	 * \brief Computes extended primitive variable vector
	 * @param V the primitive vector
	 * @return The vector (vapour mass concentration, pressure, mixture velocity, temperature, vapour volume fraction, relative velocity, vapour density, liquid density mixture enthalpy) */
	Vector computeExtendedPrimState(double *V);

	/** \fn relative_velocity
	 * \brief Computes the relative velocity between the two phases
	 * @param c_v is the vapour mass concentration: alpha_v rho_v/(alpha_v rho_v + alpha_l rho_l)
	 * @param mean_velocity is the mixture velocity (alpha_v rho_v u_v+ alpha_l rho_l u_l)/(alpha_v rho_v + alpha_l rho_l)
	 * @param rho_m is the mixture density : alpha_v rho_v + alpha_l rho_l
	 * @return The relative velocity vector u_v-u_l */

	Vector relative_velocity(double c_v, Vector mean_velocity, double rho_m)
	{
		Vector result(mean_velocity.getNumberOfRows());
		for(int i=0;i<mean_velocity.getNumberOfRows();i++)
			result(i)=0;
		return result;
	}
	/** \fn getMixtureTemperature
	 * \brief Computes the temperature of the two phase mixture from its pressure, mixture density and Gas concentration
	 * @param c_v is the vapour mass concentration: alpha_v rho_v/(alpha_v rho_v + alpha_l rho_l)
	 * @param pressureure is the mixture pressure
	 * @param rho_m is the mixture density : alpha_v rho_v + alpha_l rho_l
	 * @return The mixture temperature*/
	double getMixtureTemperature(double c_v, double rho_m, double pressure);

	/** \fn getMixturePressure
	 * \brief Computes the pressure of the two phase mixture from its temperature, mixture density and Gas concentration
	 * @param c_v is the vapour mass concentration: alpha_v rho_v/(alpha_v rho_v + alpha_l rho_l)
	 * @param temperature is the mixture temperature
	 * @param rho_m is the mixture density : alpha_v rho_v + alpha_l rho_l
	 * @return The mixture pressure */
	double getMixturePressure(double c_v, double rho_m, double temperature);

	/** \fn getMixturePressureAndTemperature
	 * \brief Computes the pressure and temperature of the two phase mixture from its mixture density, Gas concentration and internal energy
	 * @param c_v is the vapour mass concentration: alpha_v rho_v/(alpha_v rho_v + alpha_l rho_l)
	 * @param rho_m is the mixture density : alpha_v rho_v + alpha_l rho_l
	 * @param rhom_em is the mixture internal energy
	 * @return P and T */
	void getMixturePressureAndTemperature(double c_v, double rho_m, double rhom_em, double& P, double& T);

	/** \fn getMixturePressureDerivatives
	 * \brief Computes the pressure of the two phase mixture from its temperature, mixture density and Gas concentration
	 * @param m_v is the vapour partial density: alpha_v rho_v)
	 * @param m_l is the vapour partial density: alpha_l rho_l)
	 * @param pressure is the mixture pressure
	 * @param temperature is the mixture temperature
	*/
	void getMixturePressureDerivatives(double m_v, double m_l, double pression, double Tm);
};
#endif /* DRIFTMODEL_HXX_*/



