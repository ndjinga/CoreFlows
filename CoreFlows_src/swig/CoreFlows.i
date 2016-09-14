%module CoreFlows

%include std_string.i
%include std_vector.i

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
#include "SinglePhase.hxx"

%}

%include "ProblemCoreFlows.hxx"
%include "ProblemFluid.hxx"
%include "DriftModel.hxx"
%include "FiveEqsTwoFluid.hxx"
%include "IsothermalTwoFluid.hxx"
%include "TransportEquation.hxx"
%include "DiffusionEquation.hxx"
%include "SinglePhase.hxx"

