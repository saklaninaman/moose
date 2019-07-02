//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ISOTROPICDAMAGE_H
#define ISOTROPICDAMAGE_H

#include "ColumnMajorMatrix.h"
#include "ScalarDamageBase.h"
#include "SmearedCrackSofteningBase.h"
#include "Function.h"
#include "GuaranteeConsumer.h"

class IsotropicDamage;

template <>
InputParameters validParams<IsotropicDamage>();

/**
 * IsotropicDamage computes the stress for a finite strain
 * material with smeared cracking
 */
class IsotropicDamage : public ScalarDamageBase, public GuaranteeConsumer
{
public:
  IsotropicDamage(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initQpStatefulProperties() override;

protected:

  virtual void updateQpDamageIndex() override;
  void computeDamageEvolution(Real max_principal_strain, Real cracking_strain, Real cracking_stress, Real youngs_modulus);


  /// Enum defining the crack release model
  const enum class CrackingRelease { exponential, linear } _cracking_release;

  /// Threshold at which cracking initiates if tensile stress exceeds it
  const VariableValue & _tensile_strength;

  /// Ratio of the residual stress after being fully cracked to the tensile
  /// cracking threshold stress
  const Real _residual_frac;

  /// Fracture energy
  const Real _Gf;

  //@{ Strain upon crack initiation
  MaterialProperty<Real> & _crack_initiation_strain;
  const MaterialProperty<Real> & _crack_initiation_strain_old;
  ///@}

  //@{ Flag variable to indicate if cracking has occured or not
  MaterialProperty<Real> & _crack_flag;
  const MaterialProperty<Real> & _crack_flag_old;
  ///@}

  //@{ Cracking surface
  MaterialProperty<Real> & _cracking_yield_surface;
  const MaterialProperty<Real> & _cracking_yield_surface_old;
  ///@}

  //@{ Actual cracking stress during initiation
  MaterialProperty<Real> & _actual_cracking_stress;
  const MaterialProperty<Real> & _actual_cracking_stress_old;
  ///@}

  //@{ Equivalent Strain
  MaterialProperty<Real> & _equivalent_strain;
  const MaterialProperty<Real> & _equivalent_strain_old;
  ///@}

  /// Enum defining the equivalent strain definition
  const enum class EquivalentStrainDefinition { mazar, modifiedvonmises } _equivalent_strain_definition;

  /// Name of elasticity tensor
  const std::string _elasticity_tensor_name;

  /// Current undamaged elasticity tensor
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// Current stress
  const MaterialProperty<RankTwoTensor> & _stress;

  /// Current mechanical strain
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;

};

#endif // ISOTROPICDAMAGE_H
