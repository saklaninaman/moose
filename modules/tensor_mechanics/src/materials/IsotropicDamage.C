//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IsotropicDamage.h"
#include "MooseMesh.h"
#include "ElasticityTensorTools.h"
#include "StressUpdateBase.h"
#include "Conversion.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <fstream>

registerMooseObject("TensorMechanicsApp", IsotropicDamage);

template <>
InputParameters
validParams<IsotropicDamage>()
{
  InputParameters params = validParams<ComputeMultipleInelasticStress>();
  params.addClassDescription("Compute stress using an isotropic damage model");
  MooseEnum cracking_release("exponential linear", "linear");
  params.addParam<MooseEnum>("cracking_release",
			      cracking_release,
      			     "The cracking release type.  'linear' (default) gives an linear stress release "
                             "'exponential' uses an exponential softening model, ");
  params.addRequiredCoupledVar("cracking_stress",
      			       "The stress threshold beyond which cracking occurs. Negative values prevent cracking.");
  params.addParam<Real>("residual_fraction",
      			 0.0,
     			"The fraction of the cracking stress allowed to be maintained following a crack.");
  params.addParam<Real>("fracture_energy",
			"The fracture energy of the material");
  MooseEnum damage_model("mazar modifiedvonmises", "modifiedvonmises");
  params.addParam<MooseEnum>("damage_model",
			      damage_model,
      			     "How damage evolves with strain.  'modifiedvonmises' (default) calculate strain as modified version of Vonmises"
                             "'mazar' calculates strain as per Mazar model, ");
  return params;
}

IsotropicDamage::IsotropicDamage(const InputParameters & parameters)
  : ComputeMultipleInelasticStress(parameters),
    _cracking_release(getParam<MooseEnum>("cracking_release").getEnum<CrackingRelease>()),
    _cracking_stress(coupledValue("cracking_stress")),
    _residual_frac(getParam<Real>("residual_fraction")),
    _Gf(getParam<Real>("fracture_energy")),
    _crack_damage(declareProperty<Real>(_base_name + "crack_damage")),
    _crack_damage_old(getMaterialPropertyOld<Real>(_base_name + "crack_damage")),
    _crack_initiation_strain(declareProperty<Real>(_base_name + "crack_initiation_strain")),
    _crack_initiation_strain_old(getMaterialPropertyOld<Real>(_base_name + "crack_initiation_strain")),
    _crack_flag(declareProperty<int>(_base_name+"crack_flag")),
    _crack_flag0(declareProperty<Real>(_base_name+"crack_flag0")),
    _crack_flag_old(getMaterialPropertyOld<int>(_base_name+"crack_flag")),
    _crack_flag0_old(getMaterialPropertyOld<Real>(_base_name+"crack_flag0")),
    _cracking_yield_surface(declareProperty<Real>(_base_name + "crack_yield_surface")),
    _cracking_yield_surface_old(getMaterialPropertyOld<Real>(_base_name + "crack_yield_surface")),
    _actual_cracking_stress(declareProperty<Real>("actual_cracking_stress")),
    _actual_cracking_stress_old(getMaterialPropertyOld<Real>(_base_name + "actual_cracking_stress")),
    _equivalent_strain(declareProperty<Real>(_base_name +"equivalent_strain")),
    _equivalent_strain_old(getMaterialPropertyOld<Real>(_base_name + "equivalent_strain")),
    _damage_model(getParam<MooseEnum>("damage_model").getEnum<DamageModel>())
{
}

void
IsotropicDamage::initQpStatefulProperties()
{
  ComputeMultipleInelasticStress::initQpStatefulProperties();
  _crack_damage[_qp] = 0.0;
  _crack_initiation_strain[_qp]=0.0;
  _crack_flag[_qp]=0;
  _crack_flag0[_qp]=0.0;
  _cracking_yield_surface[_qp] =_cracking_stress[_qp];
  _actual_cracking_stress[_qp]=_cracking_stress[_qp];
  _equivalent_strain[_qp]=0.0;
}

void
IsotropicDamage::initialSetup()
{
  ComputeMultipleInelasticStress::initialSetup();

  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("IsotropicDamage requires that the elasticity tensor be "
               "guaranteed isotropic");
}

void
IsotropicDamage::computeQpStress()
{
  if (!previouslyCracked())
    computeQpStressIntermediateConfiguration();
  else
  {
    _elastic_strain[_qp] = _elastic_strain_old[_qp] + _strain_increment[_qp];

    // Propagate behavior from the (now inactive) inelastic models
    _inelastic_strain[_qp] = _inelastic_strain[_qp]+_inelastic_strain_old[_qp];
    for (auto model : _models)
    {
      model->setQp(_qp);
      model->propagateQpStatefulProperties();
    }

    // Since the other inelastic models are inactive, they will not limit the time step
    _matl_timestep_limit[_qp] = std::numeric_limits<Real>::max();

    // update _local_elasticity_tensor based on cracking state in previous time step
    updateLocalElasticityTensor();
      
    // Calculate stress in intermediate configuration
    _stress[_qp] = _local_elasticity_tensor * _elastic_strain[_qp];        
    _Jacobian_mult[_qp] = _local_elasticity_tensor;
  }

  // compute crack status and adjust stress
  updateCrackingStateAndStress();
}

void
IsotropicDamage::updateLocalElasticityTensor()
{
  const Real youngs_modulus =
      ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real poissons_ratio =
      ElasticityTensorTools::getIsotropicPoissonsRatio(_elasticity_tensor[_qp]);

    _local_elasticity_tensor = _elasticity_tensor[_qp];
  
	// A yield surface less than equal to zero implies no failure in tension
  	if (_cracking_yield_surface[_qp] > 0)
  	{
      	// Update elasticity tensor based on crack status of the end of last time step
      	if (_crack_damage_old[_qp] > 0.0)
      	{
        Real stiffness_ratio_local = std::max((1.0 - _crack_damage_old[_qp]),0.0);
      	std::vector<Real> local_elastic(2);
      	local_elastic[0]=youngs_modulus*stiffness_ratio_local;
      	local_elastic[1]=poissons_ratio;
      	_local_elasticity_tensor.fillSymmetricIsotropicEandNu(local_elastic[0],local_elastic[1]);
      	}
  	}
}

void
IsotropicDamage::updateCrackingStateAndStress()
{
  const Real youngs_modulus =
      ElasticityTensorTools::getIsotropicYoungsModulus(_elasticity_tensor[_qp]);
  const Real poissons_ratio =
      ElasticityTensorTools::getIsotropicPoissonsRatio(_elasticity_tensor[_qp]);
      std::vector<Real> eigval(3, 0.0);
	_crack_flag[_qp]=_crack_flag_old[_qp];
	_crack_flag0[_qp]=_crack_flag0_old[_qp];
	_crack_initiation_strain[_qp]=_crack_initiation_strain_old[_qp];
	_actual_cracking_stress[_qp]=_actual_cracking_stress_old[_qp];

    _stress[_qp].symmetricEigenvalues(eigval);
    Real max_principal_stress=std::max(std::max(eigval[0],eigval[1]),eigval[2]);
    // Find equivalent strain  
    _elastic_strain[_qp].symmetricEigenvalues(eigval);


  switch (_damage_model)
  {
    	case DamageModel::modifiedvonmises:
	{
    	Real k=10;
   	Real eps_xx=_elastic_strain[_qp](0,0);
   	Real eps_yy=_elastic_strain[_qp](1,1);
    	Real eps_zz=_elastic_strain[_qp](2,2);
    	Real eps_xy=_elastic_strain[_qp](0,1);
    	Real eps_xz=_elastic_strain[_qp](0,2);
    	Real eps_yz=_elastic_strain[_qp](1,2);
    	Real I1=eps_xx+eps_yy+eps_zz;
    	Real J2=(eps_xx*eps_xx+eps_yy*eps_yy+eps_zz*eps_zz-eps_xx*eps_yy-eps_xx*eps_zz-eps_yy*eps_zz+3*(eps_xy*eps_xy+eps_xz*eps_xz+eps_yz*eps_yz))/3;
        _equivalent_strain[_qp] = (k-1)/(2*k*(1-2*poissons_ratio))*I1+1/(2*k)*std::sqrt(std::pow(((k-1)*I1/(1-2*poissons_ratio)),2)+12*k*J2/(std::pow((1+poissons_ratio),2))); 
        break;
	}
	case DamageModel::mazar:
	{
   	Real eps_x1=std::max(eigval[0],0.0)*std::max(eigval[0],0.0);
        Real eps_x2=std::max(eigval[1],0.0)*std::max(eigval[1],0.0);
        Real eps_x3=std::max(eigval[2],0.0)*std::max(eigval[2],0.0);
     	_equivalent_strain[_qp] = std::sqrt(eps_x1+eps_x2+eps_x3);
        break;
	}
  }
     	// When the material cracks the first time
	if (max_principal_stress >= _cracking_yield_surface[_qp] && _crack_flag[_qp] == 0)
	{
        	_crack_initiation_strain[_qp] =_equivalent_strain[_qp];
        	_stress[_qp].symmetricEigenvalues(eigval); 
        	_actual_cracking_stress[_qp] =std::max(std::max(eigval[0],eigval[1]),eigval[2]);
		computeDamageEvolution(_equivalent_strain[_qp],_crack_initiation_strain[_qp],_actual_cracking_stress[_qp],youngs_modulus);      
		_crack_flag[_qp]=1;
		_crack_flag0[_qp]=1;
	}  
     	// If the crack is propagating 
        if (_crack_flag[_qp]==1 && _equivalent_strain[_qp]>_equivalent_strain_old[_qp])
       	{
   		computeDamageEvolution(_equivalent_strain[_qp],_crack_initiation_strain[_qp],_actual_cracking_stress[_qp],youngs_modulus); 
        }
	// If the crack exist but is not propagating
 	else
	{
		_crack_damage[_qp]=_crack_damage_old[_qp];
		_cracking_yield_surface[_qp]=_cracking_yield_surface_old[_qp];
  	}
}

void IsotropicDamage::computeDamageEvolution
               (Real equivalent_strain, Real cracking_strain, Real cracking_stress, Real youngs_modulus)
{
  Real residual_stress = _residual_frac*cracking_stress;
  // Get characteristic length for element
  Real hce=_current_elem->volume();
  if (_mesh.dimension()==2)
  {
  if (_current_elem->type() == 3 || _current_elem->type() == 4)
  hce=std::sqrt(2*_current_elem->volume());
  if (_current_elem->type() == 5 || _current_elem->type() == 6 || _current_elem->type() == 7)
  hce=std::sqrt(_current_elem->volume());
  }
  if (_mesh.dimension()==3)
  {
  if (_current_elem->type() == 8 || _current_elem->type() == 9)
  hce=std::cbrt(5*_current_elem->volume());
  if (_current_elem->type() == 10 || _current_elem->type() == 11 || _current_elem->type() == 12 )
  hce=std::cbrt(_current_elem->volume());
  }

  switch (_cracking_release)
  {
    	case CrackingRelease::linear:
	{
  	// Compute fracture strain
  	const Real fracture_strain=2*_Gf/(cracking_stress*hce);
        const Real maximum_strain=cracking_strain-(cracking_stress-residual_stress)*(cracking_strain-fracture_strain)/cracking_stress;
        if (equivalent_strain<maximum_strain)
   	_crack_damage[_qp]=std::min(std::max(1-cracking_stress*(equivalent_strain-fracture_strain)/(cracking_strain-fracture_strain)/youngs_modulus/equivalent_strain,_crack_damage_old[_qp]),1.0);
        else
        _crack_damage[_qp]=std::min(std::max(1-residual_stress/youngs_modulus/equivalent_strain,_crack_damage_old[_qp]),1.0);

        _cracking_yield_surface[_qp]=std::max(std::min(cracking_stress*(equivalent_strain-fracture_strain)/(cracking_strain-fracture_strain),cracking_stress),residual_stress);

        break;
	}
	case CrackingRelease::exponential:
	{
	const Real beta=youngs_modulus*cracking_strain*hce/_Gf;
        Real maximum_strain=99999999999;        
	if (residual_stress>0)
        maximum_strain = cracking_strain-1/beta*std::log(residual_stress/cracking_stress);
        if (equivalent_strain<maximum_strain)
        _crack_damage[_qp]=std::min(std::max(1-cracking_strain/equivalent_strain*std::exp(-beta*(equivalent_strain-cracking_strain)),_crack_damage_old[_qp]),1.0);
        else
        _crack_damage[_qp]=std::min(std::max(1-residual_stress/youngs_modulus/equivalent_strain,_crack_damage_old[_qp]),1.0);

	_cracking_yield_surface[_qp]=std::max(std::min(cracking_stress*std::exp(-beta*(equivalent_strain-cracking_strain)),cracking_stress),residual_stress);
	break;
	}
  }
}

bool
IsotropicDamage::previouslyCracked()
{
    if (_crack_damage_old[_qp] > 0.0)
      return true;
  return false;
}
