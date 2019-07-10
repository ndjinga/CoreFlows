#include "StationaryDiffusionEquation.hxx"
#include "math.h"
#include <algorithm> 
#include <fstream>
#include <sstream>

using namespace std;

int StationaryDiffusionEquation::fact(int n)
{
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}

StationaryDiffusionEquation::StationaryDiffusionEquation(int dim, bool FECalculation, double lambda){
	PetscBool petscInitialized;
	PetscInitialized(&petscInitialized);
	if(!petscInitialized)
		PetscInitialize(NULL,NULL,0,0);

    if(lambda < 0.)
    {
        std::cout<<"conductivity="<<lambda<<endl;
        throw CdmathException("Error : conductivity parameter lambda cannot  be negative");
    }
    if(dim==0)
    {
        std::cout<<"space dimension="<<dim<<endl;
        throw CdmathException("Error : cparameter dim cannot  be zero");
    }

    _FECalculation=FECalculation;
    
	_Ndim=dim;
	_nVar=1;//scalar prolem
	_dt_src=0;
	_diffusionMatrixSet=false;
	_initializedMemory=false;

    //Mesh data
    _neibMaxNbCells=0;    
    _neibMaxNbNodes=0;    
    _meshSet=false;
    _neibMaxNbNodes=0;
    _boundaryNodeIds=std::vector< int >(0,0);
    _NboundaryNodes=0;
    _NinteriorNodes=0;
    
    //Linear solver data
	_precision=1.e-6;
	_precision_Newton=_precision;
	_MaxIterLinearSolver=0;//During several newton iterations, stores the max petssc interations
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
	_fileName = "StationaryDiffusionProblem";
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
    {
		cout<<"Initialisation of the computation of the temperature diffusion in a solid using ";
        if(!_FECalculation)
            cout<< "Finite volumes method"<<endl;
        else
            cout<< "Finite elements method"<<endl;
    }
    
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
        _heatPowerFieldSet=true;
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
        _fluidTemperatureFieldSet=true;
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
	VecDuplicate(_Tk, &_b);//RHS of the linear system

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
	if(!_diffusionMatrixSet)//The diffusion matrix is computed once and for all time steps
		computeDiffusionMatrix(stop);

	_dt_src=computeRHS(stop);
	stop=false;
	return _dt_src;
}

Vector StationaryDiffusionEquation::gradientNodal(Matrix M, vector< double > values){
    vector< Matrix > matrices(_Ndim);
    
    for (int idim=0; idim<_Ndim;idim++){
        matrices[idim]=M.deepCopy();
        for (int jdim=0; jdim<_Ndim+1;jdim++)
			matrices[idim](jdim,idim) = values[jdim] ;
    }

	Vector result(_Ndim);
    for (int idim=0; idim<_Ndim;idim++)
        result[idim] = matrices[idim].determinant();

	return result;    
}

double StationaryDiffusionEquation::computeDiffusionMatrix(bool & stop)
{
    double result;
    
    if(_FECalculation)
        result=computeDiffusionMatrixFE(stop);
    else
        result=computeDiffusionMatrixFV(stop);

    //Contribution from the solid/fluid heat exchange with assumption of constant heat transfer coefficient
    //update value here if variable  heat transfer coefficient
    if(_heatTransfertCoeff>_precision)
        MatShift(_A,_heatTransfertCoeff);//Contribution from the liquit/solid heat transfer
        
    return  result;
}

double StationaryDiffusionEquation::computeDiffusionMatrixFE(bool & stop){
	Cell Cj;
	string nameOfGroup;
	double dn;
	MatZeroEntries(_A);
	VecZeroEntries(_b);
    
    Matrix M(_Ndim+1,_Ndim+1);//cell geometry matrix
    std::vector< Vector > GradShapeFuncs(_Ndim+1);//shape functions of cell nodes
    std::vector< int > nodeIds(_Ndim+1);//cell node Ids
    std::vector< Node > nodes(_Ndim+1);//cell nodes
    int i_int, j_int; //index of nodes j and k considered as interior nodes
    bool borderCell;
    
    std::vector< vector< double > > values(_Ndim+1,vector< double >(_Ndim+1,0));//values of shape functions on cell node
    for (int idim=0; idim<_Ndim+1;idim++)
        values[idim][idim]=1;

    /* parameters for boundary treatment */
    vector< double > valuesBorder(_Ndim+1);
    Vector GradShapeFuncBorder(_Ndim+1);
    
	for (int j=0; j<_Nmailles;j++)
    {
		Cj = _mesh.getCell(j);

        for (int idim=0; idim<_Ndim+1;idim++){
            nodeIds[idim]=Cj.getNodeId(idim);
            nodes[idim]=_mesh.getNode(nodeIds[idim]);
            for (int jdim=0; jdim<_Ndim;jdim++)
                M(idim,jdim)=nodes[idim].getPoint()[jdim];
            M(idim,_Ndim)=1;
        }
        for (int idim=0; idim<_Ndim+1;idim++)
            GradShapeFuncs[idim]=gradientNodal(M,values[idim])/fact(_Ndim);
            
        for (int idim=0; idim<_Ndim+1;idim++)
        {
            if(find(_boundaryNodeIds.begin(),_boundaryNodeIds.end(),nodeIds[idim])==_boundaryNodeIds.end()) //or for better performance nodeIds[idim]>boundaryNodes.upper_bound()
            {
                i_int=nodeIds[idim]-_NboundaryNodes;//assumes node numbering starts with interior nodes. otherwise interiorNodes.index(j)
                borderCell=false;
                for (int jdim=0; jdim<_Ndim+1;jdim++)
                {
                    if(find(_boundaryNodeIds.begin(),_boundaryNodeIds.end(),nodeIds[jdim])==_boundaryNodeIds.end()) //or for better performance nodeIds[jdim]>boundaryNodes.upper_bound()
                    {
                        j_int= nodeIds[jdim]-_NboundaryNodes;
                        MatSetValue(_A,i_int,j_int,_conductivity*(_DiffusionTensor*GradShapeFuncs[idim])*GradShapeFuncs[jdim]/Cj.getMeasure(), ADD_VALUES);
                    }
                    else if (!borderCell)
                    {
                        borderCell=true;
                        for (int kdim=0; kdim<_Ndim+1;kdim++)
                        {
                            if(nodes[kdim].isBorder())
                            {
                                nameOfGroup = nodes[kdim].getGroupName();
                                valuesBorder[kdim]=_limitField[nameOfGroup].T;
                            }
                            else
                                valuesBorder[kdim]=0;                            
                        }
                        GradShapeFuncBorder=gradientNodal(M,valuesBorder)/fact(_Ndim);
                        double coeff =-_conductivity*(_DiffusionTensor*GradShapeFuncBorder)*GradShapeFuncs[idim]/Cj.getMeasure();
                        VecSetValue(_b,i_int,coeff, ADD_VALUES);                        
                    }
                }
            }
        }            
	}
    
    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(_b);
	VecAssemblyEnd(_b);

	_diffusionMatrixSet=true;
    stop=false ;

	return INFINITY;
}

double StationaryDiffusionEquation::computeDiffusionMatrixFV(bool & stop){
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
	VecZeroEntries(_b);
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
				VecSetValue(_b,idm,inv_dxi*dn*_limitField[nameOfGroup].T, ADD_VALUES);
			}
			else {
                stop=true ;
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
	VecAssemblyBegin(_b);
	VecAssemblyEnd(_b);
    
	_diffusionMatrixSet=true;
    stop=false ;
	
    return INFINITY;
}

double StationaryDiffusionEquation::computeRHS(bool & stop)//Contribution of the PDE RHS ti the linear systemm RHS (boundary conditions do contribute to the system RHS)
{
	VecAssemblyBegin(_b);

    if(!_FECalculation)
        for (int i=0; i<_Nmailles;i++)
        {
            VecSetValue(_b,i,_heatTransfertCoeff*_fluidTemperatureField(i),ADD_VALUES);
            VecSetValue(_b,i,_heatPowerField(i)                           ,ADD_VALUES);
        }
    else
        {
            Cell Ci;
            IntTab nodesId;
            for (int i=0; i<_Nmailles;i++)
            {
                Ci=_mesh.getCell(i);
                nodesId=Ci.getNodesId();
                for (int j=0; j<nodesId.size();j++)
                    if(find(_boundaryNodeIds.begin(),_boundaryNodeIds.end(),nodesId[j])==_boundaryNodeIds.end()) //or for better performance nodeIds[idim]>boundaryNodes.upper_bound()
                    {
                        VecSetValue(_b,nodesId[j]-_NboundaryNodes,_heatTransfertCoeff*_fluidTemperatureField(nodesId[j])*Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
                        VecSetValue(_b,nodesId[j]-_NboundaryNodes,_heatPowerField(nodesId[j])                           *Ci.getMeasure()/(_Ndim+1),ADD_VALUES);
                    }
            }
        }
    
	VecAssemblyEnd(_b);

    if(_verbose or _system)
        VecView(_b,PETSC_VIEWER_STDOUT_SELF);

    stop=false ;
	return INFINITY;
}

bool StationaryDiffusionEquation::iterateNewtonStep(bool &converged)
{
	bool stop;

    //Only implicit scheme considered
    MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(  _A, MAT_FINAL_ASSEMBLY);

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
        converged = false;
        stop = true;
    }
    else{
        VecCopy(_Tk, _deltaT);//ici on a deltaT=Tk
        VecAXPY(_deltaT,  -1, _Tkm1);//On obtient deltaT=Tk-Tkm1

        _erreur_rel= 0;
        double Ti, dTi;

        VecAssemblyBegin(_Tk);
        VecAssemblyEnd(  _Tk);
        VecAssemblyBegin(_deltaT);
        VecAssemblyEnd(  _deltaT);

        if(!_FECalculation)
            for(int i=0; i<_Nmailles; i++)
            {
                VecGetValues(_deltaT, 1, &i, &dTi);
                VecGetValues(_Tk, 1, &i, &Ti);
                if(_erreur_rel < fabs(dTi/Ti))
                    _erreur_rel = fabs(dTi/Ti);
            }
        else
            for(int i=0; i<_NinteriorNodes; i++)
            {
                VecGetValues(_deltaT, 1, &i, &dTi);
                VecGetValues(_Tk, 1, &i, &Ti);
                if(_erreur_rel < fabs(dTi/Ti))
                    _erreur_rel = fabs(dTi/Ti);
            }
        stop=false;
        converged = (_erreur_rel <= _precision) ;//converged=convergence des iterations de Newton
    }

	VecCopy(_Tk, _Tkm1);

	return stop;
}

void StationaryDiffusionEquation::setMesh(const Mesh &M)
{
	if(_Ndim != M.getSpaceDimension() or _Ndim!=M.getMeshDimension()){
        cout<< "Problem : dim = "<<_Ndim<< " but mesh dim= "<<M.getMeshDimension()<<", mesh space dim= "<<M.getSpaceDimension()<<endl;
		*_runLogFile<<"StationaryDiffusionEquation::setMesh: mesh has incorrect dimension"<<endl;
		_runLogFile->close();
		throw CdmathException("StationaryDiffusionEquation::setMesh: mesh has incorrect space dimension");
	}

	_mesh=M;
	_Nmailles = _mesh.getNumberOfCells();
	_Nnodes =   _mesh.getNumberOfNodes();
    
	// find  maximum nb of neibourghs
    if(!_FECalculation)
    {
    	_VV=Field ("Temperature", CELLS, _mesh, 1);
        _neibMaxNbCells=_mesh.getMaxNbNeighbours(CELLS);
    }
    else
    {
        if(_Ndim==3 and not M.isTetrahedral())
        {
            cout<<"Dimension is "<<_Ndim<< ", mesh should be tetrahedral"<<endl;
            throw CdmathException("StationaryDiffusionEquation::setMesh: mesh has incorrect cell types");
        }
        if(_Ndim==2 and not M.isTriangular())
        {
            cout<<"Dimension is "<<_Ndim<< ", mesh should be triangular"<<endl;
            throw CdmathException("StationaryDiffusionEquation::setMesh: mesh has incorrect cell types");
        }

		_VV=Field ("Temperature", NODES, _mesh, 1);

        _neibMaxNbNodes=_mesh.getMaxNbNeighbours(NODES);
        _boundaryNodeIds=_mesh.getBoundaryNodeIds();
        _NboundaryNodes=_boundaryNodeIds.size();
        _NinteriorNodes=_Nnodes - _NboundaryNodes;
    }

	_meshSet=true;
}

void StationaryDiffusionEquation::setLinearSolver(linearSolver kspType, preconditioner pcType)
{
	//_maxPetscIts=maxIterationsPetsc;
	// set linear solver algorithm
	if (kspType==GMRES)
		_ksptype = (char*)&KSPGMRES;
	else if (kspType==BICGSTAB)
		_ksptype = (char*)&KSPBCGS;
	else {
		cout << "only 'GRMES' or 'BICGSTAB' is acceptable as a linear solver !!!" << endl;
		*_runLogFile << "only 'GRMES' or 'BICGSTAB' is acceptable as a linear solver !!!" << endl;
		_runLogFile->close();
		throw CdmathException("only 'GRMES' or 'BICGSTAB' algorithm is acceptable !!!");
	}
	// set preconditioner
	if (pcType == NONE){
		_pctype = (char*)&PCNONE;
	} else if (pcType ==LU){
		_pctype = (char*)&PCLU;
	} else if (pcType == ILU){
		_pctype = (char*)&PCILU;
	} else {
		cout << "only 'NONE' or 'LU' or 'ILU' preconditioners are acceptable !!!" << endl;
		*_runLogFile << "only 'NONE' or 'LU' or 'ILU' preconditioners are acceptable !!!" << endl;
		_runLogFile->close();
		throw CdmathException("only 'NONE' or 'LU' or 'ILU' preconditioners are acceptable !!!" );
	}
}

bool StationaryDiffusionEquation::solveStationaryProblem()
{
	if(!_initializedMemory)
	{
		_runLogFile->close();
		throw CdmathException("ProblemCoreFlows::run() call initialize() first");
	}
	bool stop=false; // Does the Problem want to stop (error) ?
	bool converged=false; // has the newton scheme converged (end) ?

	cout<< "Running test case "<< _fileName << " using ";
    if(!_FECalculation)
        cout<< "Finite volumes method"<<endl;
    else
        cout<< "Finite elements method"<<endl;

	_runLogFile->open((_fileName+".log").c_str(), ios::out | ios::trunc);;//for creation of a log file to save the history of the simulation
	*_runLogFile<< "Running test case "<< _fileName<<endl;

    computeDiffusionMatrix( stop);//For the moment the conductivity does not depend on the temperature (linear LHS)
    computeRHS(stop);//For the moment the heat power does not depend on the unknown temperature (linear RHS)
    if (stop){
        cout << "Error : failed computing right hand side, stopping calculation"<< endl;
        *_runLogFile << "Error : failed computing right hand side, stopping calculation"<< endl;
        throw CdmathException("Failed computing right hand side");
    }
    stop = iterateNewtonStep(converged);
    if (stop){
        cout << "Error : failed solving linear system, stopping calculation"<< endl;
        *_runLogFile << "Error : failed linear system, stopping calculation"<< endl;
        throw CdmathException("Failed solving linear system");
    }
    
    save();

	// Newton iteration loop for non linear problems
    /*
	while(!stop and !converged and _NEWTON_its<_maxNewtonIts)
	{
        computeDiffusionMatrix( stop);//case when the conductivity depends on the temperature (nonlinear LHS)
        computeRHS(stop);//case the heat power depends on the unknown temperature (nonlinear RHS)
        stop = iterateNewtonStep(converged);
        _NEWTON_its++;
	}
    if (stop){
        cout << "Error : failed solving Newton iteration "<<_NEWTON_its<<", stopping calculation"<< endl;
        *_runLogFile << "Error : failed solving Newton iteration "<<_NEWTON_its<<", stopping calculation"<< endl;
        throw CdmathException("Failed solving a Newton iteration");
    }
    else if(_NEWTON_its==_maxNewtonIts){
        cout << "Error : no convergence of Newton scheme. Maximum Newton iterations "<<_maxNewtonIts<<" reached, stopping calculation"<< endl;
        *_runLogFile << "Error : no convergence of Newton scheme. Maximum Newton iterations "<<_maxNewtonIts<<" reached, stopping calculation"<< endl;
        throw CdmathException("No convergence of Newton scheme");
    }
    else{
        cout << "Convergence of Newton scheme at iteration "<<_NEWTON_its<<", end of calculation"<< endl;
        *_runLogFile << "Convergence of Newton scheme at iteration "<<_NEWTON_its<<", end of calculation"<< endl;
        save();
    }
    */
    
	_runLogFile->close();
	return !stop;
}

void StationaryDiffusionEquation::save(){
	string resultFile(_path+"/StationaryDiffusionEquation");//Results

	resultFile+="_";
	resultFile+=_fileName;

	// create mesh and component info
    string suppress ="rm -rf "+resultFile+"_*";
    system(suppress.c_str());//Nettoyage des précédents calculs identiques
    
    if(_verbose or _system)
        VecView(_Tk,PETSC_VIEWER_STDOUT_SELF);

    double Ti; 
    if(!_FECalculation)
        for(int i=0; i<_Nmailles; i++)
            {
                VecGetValues(_Tk, 1, &i, &Ti);
                _VV(i)=Ti;
            }
    else
    {
        for(int i=0; i<_NinteriorNodes; i++)
        {
            VecGetValues(_Tk, 1, &i, &Ti);
            _VV(i+_NboundaryNodes)=Ti;//Assumes node numbering starts with border nodes
        }

        Node Ni;
        string nameOfGroup;
        for(int i=0; i<_NboundaryNodes; i++)
        {
            Ni=_mesh.getNode(i);
            if(Ni.isBorder())
            {
                nameOfGroup = Ni.getGroupName();
                cout<<"group name= "<<nameOfGroup<<endl;
                _VV(i)=_limitField[nameOfGroup].T;//Assumes node numbering starts with border nodes
            }
            else
            {
                std::cout<<"Node number i= "<<i <<" is not a boundary node. NboundaryNodes= "<< _NboundaryNodes <<endl;
                throw CdmathException("Error StationaryDiffusionEquation::save : unusual node numbering : should start with boundary nodes");
            }
        }
    }

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
void StationaryDiffusionEquation::terminate(){
	VecDestroy(&_Tk);
	VecDestroy(&_Tkm1);
	VecDestroy(&_deltaT);
	VecDestroy(&_b);
	MatDestroy(&_A);
}
