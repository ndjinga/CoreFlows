#include "DiffusionEquation.hxx"
#include "math.h"
#include <fstream>
#include <sstream>

using namespace std;

DiffusionEquation::DiffusionEquation(int dim,double rho,double cp, double lambda){
    if(_rho<_precision or cp<_precision)
    {
        std::cout<<"rho="<<_rho<<", cp= "<<cp<< ", precision= "<<_precision<<endl;
        throw CdmathException("Error : parameters rho and cp should be strictly positive");
    }
    if(lambda < 0.)
    {
        std::cout<<"conductivity="<<lambda<<endl;
        throw CdmathException("Error : conductivity parameter lambda cannot  be negative");
    }

	_conductivity=lambda;
	_cp=cp;
	_rho=rho;
	_diffusivity=_conductivity/(_rho*_cp);
	_Ndim=dim;
	_nVar=1;
	_dt_diffusion=0;
	_dt_src=0;
	_fluidTemperatureFieldSet=false;
	_fluidTemperature=0;
	_diffusionMatrixSet=false;
}

void DiffusionEquation::initialize()
{
	if(!_initialDataSet)
		throw CdmathException("DiffusionEquation::initialize() set initial data first");
	else
		cout<<"Initialising the diffusion of a solid temperature"<<endl;

	_DiffusionTensor=Matrix(_Ndim);
	for(int idim=0;idim<_Ndim;idim++)
		_DiffusionTensor(idim,idim)=1;
	/**************** Field creation *********************/

	if(!_heatPowerFieldSet){
		_heatPowerField=Field("Heat power",CELLS,_mesh,1);
		for(int i =0; i<_Nmailles; i++)
			_heatPowerField(i) = _heatSource;
	}
	if(!_fluidTemperatureFieldSet){
		_fluidTemperatureField=Field("Fluid temperature",CELLS,_mesh,1);
		for(int i =0; i<_Nmailles; i++)
			_fluidTemperatureField(i) = _fluidTemperature;
	}

	//creation de la matrice
	MatCreateSeqAIJ(PETSC_COMM_SELF, _Nmailles, _Nmailles, (1+_neibMaxNb), PETSC_NULL, &_A);
	VecCreate(PETSC_COMM_SELF, &_Tk);
	VecSetSizes(_Tk,PETSC_DECIDE,_Nmailles);
	VecSetFromOptions(_Tk);
	VecDuplicate(_Tk, &_Tn);
	VecDuplicate(_Tk, &_Tkm1);
	VecDuplicate(_Tk, &_deltaT);
	VecDuplicate(_Tk, &_b);//RHS of the linear system: _b=Tn/dt + _b0 + puisance volumique + couplage thermique avec le fluide
	VecDuplicate(_Tk, &_b0);//part of the RHS that comes from the boundary conditions

	for(int i =0; i<_Nmailles;i++)
		VecSetValue(_Tn,i,_VV(i), INSERT_VALUES);

	//Linear solver
	KSPCreate(PETSC_COMM_SELF, &_ksp);
	KSPSetType(_ksp, _ksptype);
	// if(_ksptype == KSPGMRES) KSPGMRESSetRestart(_ksp,10000);
	KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);
	KSPGetPC(_ksp, &_pc);
	PCSetType(_pc, _pctype);

	_initializedMemory=true;
	save();//save initial data
}

double DiffusionEquation::computeTimeStep(bool & stop){
	if(!_diffusionMatrixSet)//The diffusion matrix is computed once and for all time steps
    {
		_dt_diffusion=computeDiffusionMatrix();
        //Contribution from the solid/fluid heat exchange
        if(_timeScheme == Implicit and _heatTransfertCoeff/(_rho*_cp)>_precision)
        {   
            MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
            MatShift(_A,_heatTransfertCoeff/(_rho*_cp));
        }
    }

	_dt_src=computeRHS();

	stop=false;
	return min(_dt_diffusion,_dt_src);
}
double DiffusionEquation::computeDiffusionMatrix(){
	long nbFaces = _mesh.getNumberOfFaces();
	Face Fj;
	Cell Cell1,Cell2;
	string nameOfGroup;
	double inv_dxi, inv_dxj, inv_dxij;
	Vector normale(_Ndim);
	double dn;
	PetscInt idm, idn;
	IntTab idCells;
	MatZeroEntries(_A);
	VecZeroEntries(_b0);
	for (int j=0; j<nbFaces;j++){
		Fj = _mesh.getFace(j);

		// compute the normal vector corresponding to face j : from idCells[0] to idCells[1]
		idCells = Fj.getCellsId();
		Cell1 = _mesh.getCell(idCells[0]);
		idm = idCells[0];
		if (_Ndim >1){
			for(int l=0; l<Cell1.getNumberOfFaces(); l++){
				if (j == Cell1.getFacesId()[l]){
					for (int idim = 0; idim < _Ndim; ++idim)
						normale[idim] = Cell1.getNormalVector(l,idim);
					break;
				}
			}
		}else{ // _Ndim = 1 : assume that this is normal mesh : the face index increases in positive direction
			if (Fj.getNumberOfCells()<2) {
				if (j==0)
					normale[0] = -1;
				else if (j==nbFaces-1)
					normale[0] = 1;
				else
					throw CdmathException("DiffusionEquation::ComputeTimeStep(): computation of normal vector failed");
			} else if(Fj.getNumberOfCells()==2){
				if (idCells[0] < idCells[1])
					normale[0] = 1;
				else
					normale[0] = -1;
			}
		}
		//Compute velocity at the face Fj
		dn=_diffusivity*(_DiffusionTensor*normale)*normale;
		if(fabs(dn)>_maxvp)
			_maxvp=fabs(dn);

		// compute 1/dxi = volume of Ci/area of Fj
		if (_Ndim > 1)
			inv_dxi = Fj.getMeasure()/Cell1.getMeasure();
		else
			inv_dxi = 1/Cell1.getMeasure();

		// If Fj is on the boundary
		if (Fj.getNumberOfCells()==1) {
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "face numero " << j << " cellule frontiere " << idCells[0] << " ; vecteur normal=(";
				for(int p=0; p<_Ndim; p++)
					cout << normale[p] << ",";
				cout << ") "<<endl;
			}
			nameOfGroup = Fj.getGroupName();

			if (_limitField[nameOfGroup].bcType==Neumann){//Nothing to do
			}
			else if(_limitField[nameOfGroup].bcType==Dirichlet){
				MatSetValue(_A,idm,idm,inv_dxi*dn, ADD_VALUES);
				VecSetValue(_b0,idm,inv_dxi*dn*_limitField[nameOfGroup].T, ADD_VALUES);
			}
			else {
				cout<<"Boundary condition not treated for boundary named "<<nameOfGroup<<endl;
				cout<<"Accepted boundary condition are Neumann and Dirichlet "<<nameOfGroup<<endl;
				throw CdmathException("Unknown boundary condition");
			}

			// if Fj is inside the domain
		} else 	if (Fj.getNumberOfCells()==2 ){
			if(_verbose && _nbTimeStep%_freqSave ==0)
			{
				cout << "face numero " << j << " cellule gauche " << idCells[0] << " cellule droite " << idCells[1];
				cout << " ; vecteur normal=(";
				for(int p=0; p<_Ndim; p++)
					cout << normale[p] << ",";
				cout << ") "<<endl;
			}
			Cell2 = _mesh.getCell(idCells[1]);
			idn = idCells[1];
			if (_Ndim > 1)
				inv_dxj = Fj.getMeasure()/Cell2.getMeasure();
			else
				inv_dxj = 1/Cell2.getMeasure();
			inv_dxij=2/(1/inv_dxi+1/inv_dxj);

			MatSetValue(_A,idm,idm,inv_dxij*dn, ADD_VALUES);
			MatSetValue(_A,idm,idn,-inv_dxij*dn, ADD_VALUES);
			MatSetValue(_A,idn,idn,inv_dxij*dn, ADD_VALUES);
			MatSetValue(_A,idn,idm,-inv_dxij*dn, ADD_VALUES);
		}
		else
			throw CdmathException("DiffusionEquation::ComputeTimeStep(): incompatible number of cells around a face");
	}
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b0);
	VecAssemblyEnd(_b0);
	_diffusionMatrixSet=true;
	cout<<"_maxvp= "<<_maxvp<< " cfl= "<<_cfl<<" minl= "<<_minl<<endl;
	if(fabs(_maxvp)<_precision)
		throw CdmathException("DiffusionEquation::computeDiffusionMatrix(): maximum eigenvalue for time step is zero");
	else
		return _cfl*_minl*_minl/_maxvp;
}
double DiffusionEquation::computeRHS(){
	VecCopy(_b0,_b);

	VecAssemblyBegin(_b);      
    double Ti;  
	for (int i=0; i<_Nmailles;i++)
    {
		VecSetValue(_b,i,_heatPowerField(i)/(_rho*_cp),ADD_VALUES);//Contribution of the volumic heat power
        //Contribution due to fluid/solide heat exchange
        if(_timeScheme == Explicit)
        {
            VecGetValues(_Tn, 1, &i, &Ti);
            VecSetValue(_b,i,_heatTransfertCoeff/(_rho*_cp)*(_fluidTemperatureField(i)-Ti),ADD_VALUES);
        }
        else//Implicit scheme    
            VecSetValue(_b,i,_heatTransfertCoeff/(_rho*_cp)* _fluidTemperatureField(i)    ,ADD_VALUES);
	}
	VecAssemblyEnd(_b);

    if(_heatTransfertCoeff>_precision)
        return _rho*_cp/_heatTransfertCoeff;
    else
        return INFINITY;
}

bool DiffusionEquation::initTimeStep(double dt){
    VecAXPY(_b, 1/_dt, _Tn);
        
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

	if(_timeScheme == Implicit)
		MatShift(_A,1/_dt);

	_dt = dt;

	if(_verbose && _nbTimeStep%_freqSave ==0)
		MatView(_A,PETSC_VIEWER_STDOUT_SELF);

	return _dt>0;
}

void DiffusionEquation::abortTimeStep(){
    //Remove contribution od dt to the RHS
	VecAXPY(_b,  -1/_dt, _Tn);
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);
    //Remove contribution od dt to the matrix
	if(_timeScheme == Implicit)
		MatShift(_A,-1/_dt);
	_dt = 0;
}

bool DiffusionEquation::iterateTimeStep(bool &converged)
{
	bool stop=false;

	if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
		_maxvp=0;
		computeTimeStep(stop);
	}
	if(stop){
		converged=false;
		return false;
	}

	if(_timeScheme == Explicit)
	{
		MatMult(_A, _Tn, _Tk);
		VecAXPY(_Tk, -1, _b);
		VecScale(_Tk, -_dt);

		converged = true;
	}
	else
	{
		MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

#if PETSC_VERSION_GREATER_3_5
		KSPSetOperators(_ksp, _A, _A);
#else
		KSPSetOperators(_ksp, _A, _A,SAME_NONZERO_PATTERN);
#endif

		if(_conditionNumber)
			KSPSetComputeEigenvalues(_ksp,PETSC_TRUE);
		KSPSolve(_ksp, _b, _Tk);

		KSPGetIterationNumber(_ksp, &_PetscIts);
		if( _MaxIterLinearSolver < _PetscIts)
			_MaxIterLinearSolver = _PetscIts;
		if(_PetscIts>=_maxPetscIts)
		{
			cout<<"Systeme lineaire : pas de convergence de Petsc. Itérations maximales "<<_maxPetscIts<<" atteintes"<<endl;
			converged=false;
			return false;
		}
		else{
			VecCopy(_Tk, _deltaT);//ici on a deltaT=Tk
			VecAXPY(_deltaT,  -1, _Tkm1);//On obtient deltaT=Tk-Tkm1
			_erreur_rel= 0;
			double Ti, dTi;

			for(int i=0; i<_Nmailles; i++)
			{
				VecGetValues(_deltaT, 1, &i, &dTi);
				VecGetValues(_Tk, 1, &i, &Ti);
				if(_erreur_rel < fabs(dTi/Ti))
					_erreur_rel = fabs(dTi/Ti);
			}
			converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
		}
	}

	VecCopy(_Tk, _Tkm1);

	return true;
}
void DiffusionEquation::validateTimeStep()
{
	VecCopy(_Tk, _deltaT);//ici Tk=Tnp1 donc on a deltaT=Tnp1
	VecAXPY(_deltaT,  -1, _Tn);//On obtient deltaT=Tnp1-Tn

	_erreur_rel= 0;
	double Ti, dTi;

	for(int i=0; i<_Nmailles; i++)
	{
		VecGetValues(_deltaT, 1, &i, &dTi);
		VecGetValues(_Tk, 1, &i, &Ti);
		if(_erreur_rel < fabs(dTi/Ti))
			_erreur_rel = fabs(dTi/Ti);
	}
	_isStationary =(_erreur_rel <_precision);

	VecCopy(_Tk, _Tn);
	VecCopy(_Tk, _Tkm1);

	if(_verbose && _nbTimeStep%_freqSave ==0)
		cout <<"Valeur propre locale max: " << _maxvp << endl;

    //Remove the contribution from dt to prepare for new initTimeStep. The contribution from diffusion is not recomputed
	if(_timeScheme == Implicit)
		MatShift(_A,-1/_dt);
    //No need to remove the contribution to the right hand side since it is recomputed from scratch at each time step
    
	_time+=_dt;
	_nbTimeStep++;
	if (_nbTimeStep%_freqSave ==0 || _isStationary || _time>=_timeMax || _nbTimeStep>=_maxNbOfTimeStep)
        save();
}

void DiffusionEquation::save(){
	string resultFile(_path+"/DiffusionEquation");//Results

	resultFile+="_";
	resultFile+=_fileName;

    //On remplit le champ
    double Ti;
   	for(int i =0; i<_Nmailles;i++)
	{
        VecGetValues(_Tn, 1, &i, &Ti);
		_VV(i)=Ti;
    }
	_VV.setTime(_time,_nbTimeStep);

	// create mesh and component info
	if (_nbTimeStep ==0){
		string suppress ="rm -rf "+resultFile+"_*";
		system(suppress.c_str());//Nettoyage des précédents calculs identiques
        
        _VV.setInfoOnComponent(0,"Temperature_(K)");
		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(resultFile);
			break;
		case MED :
			_VV.writeMED(resultFile);
			break;
		case CSV :
			_VV.writeCSV(resultFile);
			break;
		}
	}
	else{	// do not create mesh
		switch(_saveFormat)
		{
		case VTK :
			_VV.writeVTK(resultFile,false);
			break;
		case MED :
			_VV.writeMED(resultFile,false);
			break;
		case CSV :
			_VV.writeCSV(resultFile);
			break;
		}
	}
    
    if(_isStationary)
	{
        resultFile+="_Stat";
        switch(_saveFormat)
        {
        case VTK :
            _VV.writeVTK(resultFile);
            break;
        case MED :
            _VV.writeMED(resultFile);
            break;
        case CSV :
            _VV.writeCSV(resultFile);
            break;
        }
    }
}

void DiffusionEquation::terminate(){
	VecDestroy(&_Tn);
	VecDestroy(&_Tk);
	VecDestroy(&_Tkm1);
	VecDestroy(&_deltaT);
	VecDestroy(&_b0);
	VecDestroy(&_b);
	MatDestroy(&_A);
}

