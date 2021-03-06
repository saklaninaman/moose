//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledVarNeumannBC.h"

registerMooseObject("MooseApp", CoupledVarNeumannBC);

defineLegacyParams(CoupledVarNeumannBC);

InputParameters
CoupledVarNeumannBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredCoupledVar("v", "Coupled variable setting the gradient on the boundary.");
  params.addParam<Real>(
      "coef", 1.0, "Coefficent ($\\sigma$) multiplier for the coupled force term.");
  params.addClassDescription("Imposes the integrated boundary condition "
                             "$\\frac{\\partial u}{\\partial n}=v$, "
                             "where $v$ is a variable.");
  return params;
}

CoupledVarNeumannBC::CoupledVarNeumannBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _coupled_var(coupledValue("v")), _coef(getParam<Real>("coef"))
{
}

Real
CoupledVarNeumannBC::computeQpResidual()
{
  return -_coef * _test[_i][_qp] * _coupled_var[_qp];
}
