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
#include "ComputeMultipleInelasticStress.h"
#include "SmearedCrackSofteningBase.h"
#include "Function.h"

class IsotropicDamage;

template <>
InputParameters validParams<IsotropicDamage>();

/**
 * IsotropicDamage computes the stress for a finite strain
 * material with smeared cracking
 */
class IsotropicDamage : public ComputeMultipleInelasticStress
{
public:
  IsotropicDamage(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

protected:
  /**
   * Update the local elasticity tensor (_local_elasticity_tensor)
   * due to the effects of cracking.
   */
  void updateLocalElasticityTensor();

  /**
   * Update all cracking-related state variables and the stress
   * tensor due to cracking in all directions.
   */
  virtual void updateCrackingStateAndStress();


  void computeDamageEvolution(Real max_principal_strain, Real cracking_strain, Real cracking_stress, Real youngs_modulus);

  /**
   * Check to see whether there was cracking in the previous
   * time step.
   * @return true if cracked, false if not cracked
   */
  bool previouslyCracked();

  /// Enum defining the crack release model
  const enum class CrackingRelease { exponential, linear } _cracking_release;

  /// Threshold at which cracking initiates if tensile stress exceeds it
  const VariableValue & _cracking_stress;

  /// Ratio of the residual stress after being fully cracked to the tensile
  /// cracking threshold stress
  const Real _residual_frac;

  const Real _Gf;

  //@{ Damage (goes from 0 to 1) in crack directions
  MaterialProperty<Real> & _crack_damage;
  const MaterialProperty<Real> & _crack_damage_old;
  ///@}

  //@{ Strain in direction of crack upon crack initiation
  MaterialProperty<Real> & _crack_initiation_strain;
  ///@}
  const MaterialProperty<Real> & _crack_initiation_strain_old;
  MaterialProperty<int> & _crack_flag;
  MaterialProperty<Real> & _crack_flag0;
  const MaterialProperty<int> & _crack_flag_old;
  const MaterialProperty<Real> & _crack_flag0_old;
  MaterialProperty<Real> & _cracking_yield_surface;
  const MaterialProperty<Real> & _cracking_yield_surface_old;
  MaterialProperty<Real> & _actual_cracking_stress;
  const MaterialProperty<Real> & _actual_cracking_stress_old;
  MaterialProperty<Real> & _equivalent_strain;
  const MaterialProperty<Real> & _equivalent_strain_old;
  //@{ Variables used by multiple methods within the calculation for a single material point
  RankFourTensor _local_elasticity_tensor;
  ///@}

  /// Enum defining the damage model
  const enum class DamageModel { mazar, modifiedvonmises } _damage_model;

};

#endif // ISOTROPICDAMAGE_H
