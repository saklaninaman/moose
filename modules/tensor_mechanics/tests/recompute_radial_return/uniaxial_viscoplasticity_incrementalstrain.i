# This is a test of the HyperbolicViscoplasticityStressUpdate model
# using the small strain formulation. The material is a visco-plastic material
# i.e. a time-dependent linear strain hardening plasticity model.
# A similar problem was run in Abaqus with exactly the same result, although the element
# used in the Abaqus simulation was a CAX4 element.  Neverthless, due to the boundary conditions
# and load, the MOOSE and Abaqus result are the same.

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  order = CONSTANT
  family = MONOMIAL
[]

[Mesh]
  file = 1x1x1cube.e
[]

[AuxVariables]
  [./stress_yy]
  [../]
  [./plastic_strain_xx]
  [../]
  [./plastic_strain_yy]
  [../]
  [./plastic_strain_zz]
  [../]
[]

[Functions]
  [./top_pull]
    type = ParsedFunction
    value = t/100
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = INCREMENTAL
    add_variables = true
  [../]
[]

[AuxKernels]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]

  [./plastic_strain_yy]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = plastic_strain_yy
    index_i = 1
    index_j = 1
  [../]

  [./plastic_strain_xx]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = plastic_strain_xx
    index_i = 0
    index_j = 0
  [../]

  [./plastic_strain_zz]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = plastic_strain_zz
    index_i = 2
    index_j = 2
  [../]
 []


[BCs]
  [./y_pull_function]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 5
    function = top_pull
  [../]

  [./x_bot]
    type = DirichletBC
    variable = disp_x
    boundary = 4
    value = 0.0
  [../]

  [./y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = 3
    value = 0.0
  [../]

  [./z_bot]
    type = DirichletBC
    variable = disp_z
    boundary = 2
    value = 0.0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1000.0
    poissons_ratio = 0.3
  [../]
  [./viscoplasticity]
    type = HyperbolicViscoplasticityStressUpdate
    yield_stress = 10.0
    hardening_constant = 100.0
    c_alpha = 0.2418e-6
    c_beta = 0.1135
    relative_tolerance = 1e-25
    absolute_tolerance = 1e-5
  [../]
  [./radial_return_stress]
    type = ComputeReturnMappingStress
    return_mapping_models = 'viscoplasticity'
  [../]
[]

[Executioner]
  type = Transient

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  line_search = none


  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'

  l_max_its = 100
  nl_max_its = 100
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-10
  l_tol = 1e-9

  start_time = 0.0
  num_steps = 30

  dt = 1.0
[]

[Outputs]
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
