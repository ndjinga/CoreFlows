CoreFlows example scripts
=========================

- Transport equation examples
    - [Problem Transport 1D Heated Channel (C)](../examples/TransportEquation_1DHeatedChannel.cxx) 
    - [Problem Transport 1D Heated Channel (python)](../examples/Python/TransportEquation_1DHeatedChannel.py)

- Diffusion equation examples
    - [Problem Diffusion 1D Heated Rod (C)](../examples/DiffusionEquation_1DHeatedRod.cxx)
    - [Problem Diffusion 1D Heated Rod (Python)](../examples/Python/DiffusionEquation_1DHeatedRod.py)

- Coupling of a transport and a diffusion equations
    - [Coupled Problem Transport Diffusion 1D Heated Channel (C)](../examples/CoupledTransportDiffusionEquations_1DHeatedChannel.cxx)

- Single phase examples
    - [Single Phase Problem 1D Heated Channel (C)](../examples/SinglePhase_1DHeatedChannel.cxx)
    - [Single Phase Problem 1D Heated Channel (Python)](../examples/Python/SinglePhase_1DHeatedChannel.py)
    - [Single Phase Problem 1D RiemannProblem (C)](../examples/SinglePhase_1DRiemannProblem.cxx)
    - [Single Phase Problem 1D RiemannProblem (Python)](../examples/Python/SinglePhase_1DRiemannProblem.py)
    - [Single Phase Problem 1D depressurisation (C)](../examples/SinglePhase_1DDepressurisation.cxx)
    - [Single Phase Problem 1D depressurisation (Python)](../examples/Python/SinglePhase_1DDepressurisation.py
)
    - [Single Phase Problem 1D water hammer (Python)](../examples/Python/SinglePhase_1DWaterHammer.py)
    - [Single Phase Problem 1D porosity jump (C)](../examples/SinglePhase_1DPorosityJump.cxx)

    2D Tests meshes | 
    --- | ---
    <img src="Figures/BoiteStruct.png" width="330"/> | <img src="Figures/BoiteNStruct.png" width="330"/> 
    --- | ---
    Structured mesh | Unstructured mesh

    * 2D Tests scripts
    - [Single Phase Problem 2D Lid driven cavity (C)](../examples/SinglePhase_2DLidDrivenCavity.cxx)
    - [Single Phase Problem 2D Lid driven cavity (Python)](../examples/Python/SinglePhase_2DLidDrivenCavity.py)
    - [Single Phase Problem 2D Lid driven cavity on an unstructured mesh (C)](../examples/SinglePhase_2DLidDrivenCavity_unstructured.cxx)
    - [Single Phase Problem 2D Lid driven cavity on an unstructured mesh (Python)](../examples/Python/SinglePhase_2DLidDrivenCavity_unstructured.py)

    * 2D Tests scripts results
    ![2D Lid driven cavity test with upwind scheme on a structured mesh](Figures/Simulations/DrivenCavity/DrivenCavityStructuredUpwind.png)
    ![2D Lid driven cavity test with upwind scheme on an unstructured mesh](Figures/Simulations/DrivenCavity/DrivenCavityUnstructuredUpwind.png)
    ![2D Lid driven cavity test with lowMachDellacherieOmnes scheme on a structured mesh](Figures/Simulations/DrivenCavity/DrivenCavityStructuredDellacherieOmnes.png)
    ![2D Lid driven cavity test with lowMachDellacherieOmnes scheme on an unstructured mesh](Figures/Simulations/DrivenCavity/DrivenCavityUnstructuredDellacherieOmnes.png)
    ![2D Lid driven cavity test with staggered scheme on a structured mesh](Figures/Simulations/DrivenCavity/DrivenCavityStructuredStaggered.png)
    ![2D Lid driven cavity test with staggered scheme on an unstructured mesh](Figures/Simulations/DrivenCavity/DrivenCavityUnstructuredStaggered.png)

    - [Single Phase Problem 2D Spherical wave on an unstructured mesh (C)](../examples/SinglePhase_2DSphericalExplosion_unstructured.cxx)
    - [Single Phase Problem 2D Spherical wave on an unstructured mesh (Python)](../examples/Python/SinglePhase_2DSphericalExplosion_unstructured.py)
    - [Single Phase Problem 2D Heated flow in inclined channel (Python)](../examples/Python/SinglePhase_2DHeatedChannelInclined.py)
    - [Single Phase Problem 2D Wall heated Channel with cross section change (C)](../examples/SinglePhase_2DWallHeatedChannel_ChangeSect.cxx)
    - [Single Phase Problem Heated wire with 2 Branches (C)](../examples/SinglePhase_HeatedWire_2Branches.cxx)
    
    ![Power field of the heated wire problem](Figures/2BranchesHeatedChannelPower.png)
    ![Stationary temperature field of the heated wire problem](Figures/2BranchesHeatedChannelTemperature.png)

    - [Single Phase Problem Heated wire with 2 Branches (Python)](../examples/Python/SinglePhase_2BranchesHeatedChannels.py)
    - [Single Phase Problem 2D Thermal conduction (Python)](../examples/Python/SinglePhase_2DThermalDiffusion.py)

    ![Initial data and mesh for the 2D thermal diffusion problem](Figures/2DRiemann_t0.png)
    ![Stationary solution for the 2D thermal diffusion problem](Figures/2DRiemannStaggered.png)
    
    - [Single Phase 2D tank drainage problem (Python)](../examples/Python/SinglePhase_2DVidangeReservoir.py)
    - [Single Phase Problem 3D Heat driven cavity (C)](../examples/SinglePhase_3DHeatDrivenCavity.cxx)
    - [Single Phase Problem 3D Heat driven cavity (Python)](../examples/Python/SinglePhase_3DHeatDrivenCavity.py)

- " The Drift Model"
	+ \ref Drift1DRP
	+ \ref Drift1DRPpy
	+ \ref Drift1DDepr
	+ \ref Drift1DDeprPy
	+ \ref Drift1DPressLoss
	+ \ref Drift1DPressLossPy
	+ \ref Drift1DPorJump
	+ \ref Drift1DPorJumpPy
	+ \ref Drift1DBoilChann
	+ \ref Drift1DBoilChannPy
	+ \ref Drift1DBoilAssPy
	+ \ref Drift1DVidResPy
	+ \ref Drift2DBoilC
	+ \ref Drift2DBoilPy
	+ \ref Drift2DInclGravC
	+ \ref Drift2DInclGravPy
	+ \ref Drift2DInclGravTriPy
	+ \ref Drift2DInclGravBarsC
	+ \ref Drift2DInclGravBarsPy
	+ \ref Drift2DPressLossPy
	+ \ref Drift2DPorJumpPy
	+ \ref Drift2DBoilBarVertPy
	+ \ref Drift2DBoilBarIncPy
	+ \ref Drift2DVidResPy
	+ \ref Drift2DVidResUnsPy
	+ \ref DriftWirePy 
	+ \ref Drift3DBoilBarPy

- " The Isothermal TwoFluid model "
	+ \ref IsoRP
 	+ \ref IsoSed
 	+ \ref IsoSedPy
	+ \ref IsoVidPy
	+ \ref IsoDepr
 	+ \ref IsoSed2D
	+ \ref IsoVid2D
	+ \ref IsoVid2DPy

- " The five equations model " 
	+ \ref RiemannPrblm5Eq
	+ \ref BoilingChannel5Eq
	+ \ref BoilingChannel5EqPy 
	+ \ref BoilingAssembly5EqPy 
	+ \ref FiveEqVidangeReservoirPy
	+ \ref FiveEqDepr
	+ \ref FiveEqSed2DInclined
	+ \ref FiveEqSed2DInclinedPy 
	+ \ref FiveEqBoilingChannel2D
	+ \ref FiveEqBoilingChannel2DPy 
	+ \ref FiveEqVidangeReservoir2DPy

