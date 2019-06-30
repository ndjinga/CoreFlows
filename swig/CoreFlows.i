%module CoreFlows

%include std_string.i
%include std_vector.i

%apply bool& INOUT {bool &stop}

namespace std {
 %template(VectorDouble) vector<double>;
};

%{

#include "DriftModel.hxx"
#include "FiveEqsTwoFluid.hxx"
#include "IsothermalTwoFluid.hxx"
#include "ProblemFluid.hxx"
#include "ProblemCoreFlows.hxx"
#include "TransportEquation.hxx"
#include "DiffusionEquation.hxx"
#include "StationaryDiffusionEquation.hxx"
#include "SinglePhase.hxx"
#include "Fluide.h"

%}

%include "ProblemCoreFlows.hxx"
%include "ProblemFluid.hxx"
%include "DriftModel.hxx"
%include "FiveEqsTwoFluid.hxx"
%include "IsothermalTwoFluid.hxx"
%include "TransportEquation.hxx"
%include "DiffusionEquation.hxx"
%include "StationaryDiffusionEquation.hxx"
%include "SinglePhase.hxx"
%include "Fluide.h"

