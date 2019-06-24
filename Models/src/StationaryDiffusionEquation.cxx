#include "StationaryDiffusionEquation.hxx"
#include "math.h"
#include <fstream>
#include <sstream>

using namespace std;

StationaryDiffusionEquation::StationaryDiffusionEquation(int dim, double lambda, bool FECalculation){
	PetscBool petscInitialized;
	PetscInitialized(&petscInitialized);
	if(!petscInitialized)
		PetscInitialize(NULL,NULL,0,0);

    _FECalculation=FECalculation;
    
	_Ndim=dim;
	_nVar=1;
	_dt_src=0;
	_diffusionMatrixSet=false;
    _neibMaxNbCells=0;    
    _meshSet=false;
	_initializedMemory=false;
    
    //Linear solver data
	_precision=1.e-6;
	_precision_Newton=_precision;
	_maxPetscIts=50;
	_maxNewtonIts=50;
	_NEWTON_its=0;
	int _PetscIts=0;//the number of iterations of the linear solver
	_ksptype = (char*)&KSPGMRES;
	_pctype = (char*)&PCLU;
	_conditionNumber=false;
	_erreur_rel= 0;

    //monitoring simulation
	_verbose = false;
	_system = false;
	_runLogFile=new ofstream;

	//save results
	_fileName = "myStationaryDiffusionProblem";
	char result[ PATH_MAX ];//extracting current directory
	getcwd(result, PATH_MAX );
	_path=string( result );
	_saveFormat=VTK;
    
    //heat transfer parameters
	_conductivity=lambda;
	_fluidTemperatureFieldSet=false;
	_fluidTemperature=0;
	_heatPowerFieldSet=false;
	_heatTransfertCoeff=0;
	_heatSource=0;
}

void StationaryDiffusionEquation::initialize()
{
	if(!_meshSet)
		throw CdmathException("StationaryDiffusionEquation::initialize() set initial data first");
	else
		cout<<"Initialising the diffusion of a solid temperature"<<endl;

	_DiffusionTensor=Matrix(_Ndim);
	for(int idim=0;idim<_Ndim;idim++)
		_DiffusionTensor(idim,idim)=1;
	/**************** Field creation *********************/

	if(!_heatPowerFieldSet){
        if(_FECalculation){
            _heatPowerField=Field("Heat power",NODES,_mesh,1);
            for(int i =0; i<_Nnodes; i++)
                _heatPowerField(i) = _heatSource;
        }
        else{
            _heatPowerField=Field("Heat power",CELLS,_mesh,1);
            for(int i =0; i<_Nmailles; i++)
                _heatPowerField(i) = _heatSource;
        }
    }
	if(!_fluidTemperatureFieldSet){
        if(_FECalculation){
            _fluidTemperatureField=Field("Fluid temperature",NODES,_mesh,1);
            for(int i =0; i<_Nnodes; i++)
                _fluidTemperatureField(i) = _fluidTemperature;
        }
        else{
            _fluidTemperatureField=Field("Fluid temperature",CELLS,_mesh,1);
            for(int i =0; i<_Nmailles; i++)
                _fluidTemperatureField(i) = _fluidTemperature;
        }
	}

	//creation de la matrice
    if(!_FECalculation)
        MatCreateSeqAIJ(PETSC_COMM_SELF, _Nmailles, _Nmailles, (1+_neibMaxNbCells), PETSC_NULL, &_A);
    else
        MatCreateSeqAIJ(PETSC_COMM_SELF, _NinteriorNodes, _NinteriorNodes, (1+_neibMaxNbNodes), PETSC_NULL, &_A);
	VecCreate(PETSC_COMM_SELF, &_Tk);

    if(!_FECalculation)
        VecSetSizes(_Tk,PETSC_DECIDE,_Nmailles);
    else
        VecSetSizes(_Tk,PETSC_DECIDE,_NinteriorNodes);

	VecSetFromOptions(_Tk);
	VecDuplicate(_Tk, &_Tkm1);
	VecDuplicate(_Tk, &_deltaT);
	VecDuplicate(_Tk, &_b);//RHS of the linear system: _b=Hn/dt + _b0 + puisance
	VecDuplicate(_Tk, &_b0);//part of the RHS that comes from the boundary conditions

	//Linear solver
	KSPCreate(PETSC_COMM_SELF, &_ksp);
	KSPSetType(_ksp, _ksptype);
	// if(_ksptype == KSPGMRES) KSPGMRESSetRestart(_ksp,10000);
	KSPSetTolerances(_ksp,_precision,_precision,PETSC_DEFAULT,_maxPetscIts);
	KSPGetPC(_ksp, &_pc);
	PCSetType(_pc, _pctype);

	_initializedMemory=true;

}

double StationaryDiffusionEquation::computeTimeStep(bool & stop){
	if(!_diffusionMatrixSet)
		computeDiffusionMatrix();
	_dt_src=computeRHS();

	stop=false;
	return _dt_src;
}
double StationaryDiffusionEquation::computeDiffusionMatrix(){
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
					throw CdmathException("StationaryDiffusionEquation::ComputeTimeStep(): computation of normal vector failed");
			} else if(Fj.getNumberOfCells()==2){
				if (idCells[0] < idCells[1])
					normale[0] = 1;
				else
					normale[0] = -1;
			}
		}
		//Compute velocity at the face Fj
		dn=_conductivity*(_DiffusionTensor*normale)*normale;

		// compute 1/dxi = volume of Ci/area of Fj
		if (_Ndim > 1)
			inv_dxi = Fj.getMeasure()/Cell1.getMeasure();
		else
			inv_dxi = 1/Cell1.getMeasure();

		// If Fj is on the boundary
		if (Fj.getNumberOfCells()==1) {
			if(_verbose )
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
			if(_verbose )
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
			throw CdmathException("StationaryDiffusionEquation::ComputeTimeStep(): incompatible number of cells around a face");
	}
	MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b0);
	VecAssemblyEnd(_b0);
	_diffusionMatrixSet=true;

	return INFINITY;
}
double StationaryDiffusionEquation::computeRHS(){
	VecCopy(_b0,_b);
	VecAssemblyBegin(_b);

    if(!_FECalculation)
        for (int i=0; i<_Nmailles;i++){
            VecSetValue(_b,i,_heatTransfertCoeff*(_fluidTemperatureField(i)-_VV(i)),ADD_VALUES);
            VecSetValue(_b,i,_heatPowerField(i),ADD_VALUES);
        }
    else
        for (int i=0; i<_NinteriorNodes;i++){
            VecSetValue(_b,i,_heatTransfertCoeff*(_fluidTemperatureField(i)-_VV(i)),ADD_VALUES);
            VecSetValue(_b,i,_heatPowerField(i),ADD_VALUES);
        }
    
	VecAssemblyEnd(_b);

	return INFINITY;
}

bool StationaryDiffusionEquation::iterateTimeStep(bool &converged)
{
	bool stop=false;

	if(_NEWTON_its>0){//Pas besoin de computeTimeStep à la première iteration de Newton
		computeTimeStep(stop);
	}
	if(stop){
		converged=false;
		return false;
	}
    //Only implicit scheme considered
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
            _VV(i)=Ti;
            if(_erreur_rel < fabs(dTi/Ti))
                _erreur_rel = fabs(dTi/Ti);
        }
        converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
    }

	VecCopy(_Tk, _Tkm1);

	return true;
}

void StationaryDiffusionEquation::setMesh(const Mesh &M)
{

	if(_Ndim != M.getSpaceDimension()){
		*_runLogFile<<"ProblemCoreFlows::setInitialField: mesh has incorrect space dimension"<<endl;
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::setInitialField: mesh has incorrect space dimension");
	}

	_mesh=M;
	_Nmailles = _mesh.getNumberOfCells();
	_Nnodes =   _mesh.getNumberOfNodes();
    
	// find  maximum nb of neibourghs
    if(!_FECalculation)
        _neibMaxNbCells=_mesh.getMaxNbNeighbours(CELLS);
    else
    {
        _neibMaxNbNodes=_mesh.getMaxNbNeighbours(NODES);
        _boundaryNodeIds=_mesh.getBoundaryNodeIds();
        _NboundaryNodes=_boundaryNodeIds.size();
        _NinteriorNodes=_Nnodes - _NboundaryNodes;
    }

	_meshSet=true;
}

void StationaryDiffusionEquation::save(){
	string resultFile(_path+"/StationaryDiffusionEquation");///Results

	resultFile+="_";
	resultFile+=_fileName;

	// create mesh and component info
    string suppress ="rm -rf "+resultFile+"_*";
    system(suppress.c_str());//Nettoyage des précédents calculs identiques
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
void StationaryDiffusionEquation::terminate(){
	VecDestroy(&_Tk);
	VecDestroy(&_Tkm1);
	VecDestroy(&_deltaT);
	VecDestroy(&_b0);
	VecDestroy(&_b);
	MatDestroy(&_A);
}
