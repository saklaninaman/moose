[Mesh]
  type = FileMesh
  file = PPT_ConcHum_vfpt3753_out.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./h]
    order = FIRST
    family = LAGRANGE
  [../]
  [./temp]
    order = FIRST
    family = LAGRANGE
  [../]
  [./radiation]
    order = FIRST
    family = LAGRANGE
  [../]
  [./radiationmor]
    order = FIRST
    family = LAGRANGE
  [../]
  [./time]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Isotropic_damage_index]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./MaxPrincipal_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./MaxPrincipal_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./youngs_modulus]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
[]

[AuxKernels]
  [./humidity]
    type = SolutionAux
    variable = h
    solution = soln1
  [../]
  [./tempfuncaux]
    type = SolutionAux
    variable = temp
    solution = soln2
  [../]
  [./radiationfuncaux]
    type = FunctionAux
    variable = radiation
    function = radiationfunc
    use_displaced_mesh = false
  [../]
  [./radiationmorfuncaux]
    type = FunctionAux
    variable = radiationmor
    function = radiationmorfunc
    use_displaced_mesh = false
  [../]
  [./timefuncaux]
    type = FunctionAux
    variable = time
    function = timefunc
  [../]
  [./MaxPrincipal_strain]
    type = RankTwoScalarAux
    rank_two_tensor = total_strain
    variable = MaxPrincipal_strain
    scalar_type = maxPrincipal
  [../]
  [./MaxPrincipal_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = MaxPrincipal_stress
    scalar_type = maxPrincipal
  [../]
  [./Isotropic_damage_index]
    type = MaterialRealAux
    variable = Isotropic_damage_index
    property = damage_index
  [../]
  [./youngs_modulus]
    type = MaterialRealAux
    variable = youngs_modulus
    property = youngs_modulus
    block = 1
  [../]
[]

[Postprocessors]
[]

[Functions]
  [./timefunc]
   type = ParsedFunction
   value = 't'
  [../]
  [./radiationfunc]
    type = PiecewiseLinear
    xy_data = '0 0
		2777778 0.001048616
		5555556 0.005366647
		8333333 0.01361806
		11111111 0.025615045
		13888889 0.040494315
		16666667 0.056922269
		19444444 0.073387907
		22222222 0.088519547
		25000000 0.101336688
		27777778 0.111370933'
  [../]
  [./radiationmorfunc]
    type = PiecewiseLinear
    xy_data = '0 0
345600 2.04423E-08
691200 1.06397E-07
1036800 2.79205E-07
1382400 5.53492E-07
1728000 9.40856E-07
2073600 1.451E-06
2419200 2.09224E-06
2764800 2.87183E-06
3110400 3.79612E-06
3456000 4.87066E-06
3801600 6.1003E-06
4147200 7.48923E-06
4492800 9.041E-06
4838400 1.07586E-05
5184000 1.26444E-05
5529600 1.47003E-05
5875200 1.69276E-05
6220800 1.9327E-05
6566400 2.18988E-05
6912000 2.46427E-05
7257600 2.75579E-05
7603200 3.06432E-05
7948800 3.38968E-05
8294400 3.73164E-05
8640000 4.08993E-05
8985600 4.46423E-05
9331200 4.85419E-05
9676800 5.25939E-05
10022400 5.67938E-05
10368000 6.11368E-05
10713600 6.56176E-05
11059200 7.02306E-05
11404800 7.49698E-05
11750400 7.98288E-05
12096000 8.48011E-05
12441600 8.98798E-05
12787200 9.50576E-05
13132800 0.000100327
13478400 0.000105681
13824000 0.000111112
14169600 0.00011661
14515200 0.00012217
14860800 0.000127781
15206400 0.000133437
15552000 0.000139129
15897600 0.000144848
16243200 0.000150586
16588800 0.000156336
16934400 0.000162088
17280000 0.000167836
17625600 0.00017357
17971200 0.000179283
18316800 0.000184968
18662400 0.000190617
19008000 0.000196223
19353600 0.000201778
19699200 0.000207277
20044800 0.000212712
20390400 0.000218078
20736000 0.000223368
21081600 0.000228577
21427200 0.0002337
21772800 0.000238732
22118400 0.000243669
22464000 0.000248505
22809600 0.000253238
23155200 0.000257863
23500800 0.000262378
23846400 0.00026678
24192000 0.000271066
24537600 0.000275235
24883200 0.000279284
25228800 0.000283213
25574400 0.00028702
25920000 0.000290705
26265600 0.000294267
26611200 0.000297706
26956800 0.000301022
27302400 0.000304217
27648000 0.00030729
27993600 0.000310243
28339200 0.000313077
28684800 0.000315793
29030400 0.000318394
29376000 0.000320882
29721600 0.000323257
30067200 0.000325524
30240000 0.000326616'
  [../]
[]


[Kernels]
  [./TensorMechanics]
    use_displaced_mesh = true
  [../]
[]


[BCs]
  [./zfix]
    type = DirichletBC
    variable = disp_z
    boundary = 2
    value = 0
  [../]
  [./xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0
  [../]
  [./yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0
  [../]
[]
##############################################AGGREGATE#######################################################################
[Materials]
  [./strain_aggregate]
    type = ComputeFiniteStrain
    block = 1
    eigenstrain_names = 'radiation_strain'
  [../]
  [./youngs_modulus]
    type = PiecewiseLinearInterpolationMaterial
    xy_data =  '0 72.9E9
		1388889 59.78E9
		2777778 47.96E9
		3055556 47.06E9
		3333333 46.24E9
		3611111 45.51E9
		3888889 44.84E9
		4166667 44.23E9
		4444444 43.66E9
		4722222 43.13E9
		5000000 42.64E9
		5277778 42.18E9
		5555556 41.75E9
		5833333 41.35E9
		6111111 40.96E9
		6388889 40.60E9
		6666667 40.26E9
		6944444 39.93E9
		7222222 39.62E9
		7500000 39.32E9
		7777778 39.04E9
		8055556 38.76E9
		8333333 38.50E9
		8611111 38.25E9
		8888889 38.01E9
		9166667 37.77E9
		9444444 37.55E9
		9722222 37.33E9
		10000000 37.12E9
		10277778 36.92E9
		10555556 36.72E9
		10833333 36.53E9
		11111111 36.35E9
		11388889 36.17E9
		11666667 35.99E9
		11944444 35.83E9
		12222222 35.66E9
		12500000 35.50E9
		12777778 35.35E9
		13055556 35.19E9
		13333333 35.05E9
		13611111 34.90E9
		13888889 34.76E9
		16666667 33.52E9
		19444444 32.50E9
		22222222 31.64E9
		25000000 30.91E9
		27777778 30.26E9'
    property = youngs_modulus
    variable = time
    block = 1
  [../]
  [./elasticity_tensor_aggreagte]
    type = ComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = 0.28
    args = ''
    block = 1
  [../]
  [./thermal_strain_aggregate]
    type = ComputeThermalExpansionEigenstrain
    thermal_expansion_coeff = 9e-6
    stress_free_temperature = 20.0
    temperature = temp
    incremental_form = true
    eigenstrain_name = 'thermal_strain'
    block = 1
  [../]
  [./radiation_strain_aggregate]
    type = ComputeVolumetricEigenstrain
    volumetric_materials = radiation
    eigenstrain_name = 'radiation_strain'
    args = ''
    block = 1
  [../]
  [./radiation]
    type = GenericFunctionMaterial
    prop_names = radiation
    prop_values = radiationfunc
    block = 1 
  [../]
  [./damage1]
    type = IsotropicDamage
    cracking_release = exponential
    tensile_strength = 1000E6             
    fracture_energy=0.2E3            
    residual_fraction=0.0
    equivalent_strain_definition=modifiedvonmises
    use_old_damage= 'true'
    creep_damage_parameter=0.1
    block = 1
  [../]
  [./stress1]
    type = ComputeMultipleInelasticStress
    damage_model = damage1
    inelastic_models = ''
    block = 1
  [../]
##############################################MORTAR#######################################################################
################# CREEP MODEL
  [./strain]
    type = ComputeFiniteStrain
     block = 2
     eigenstrain_names = 'radiation_strainmor'
  [../]
  [./thermal_strain_paste]
    type = ComputeThermalExpansionEigenstrain
    thermal_expansion_coeff = 11e-6
    stress_free_temperature = 20.0
    temperature = temp
    incremental_form = true
    eigenstrain_name = 'thermal_eigenstrain_paste'
    block = 2
  [../]
  [./radiation_strain_mor]
    type = ComputeVolumetricEigenstrain
    volumetric_materials = radiationmor
    eigenstrain_name = 'radiation_strainmor'
    args = ''
    block = 2
  [../]
  [./radiationmor]
    type = GenericFunctionMaterial
    prop_names = radiationmor
    prop_values = radiationmorfunc
    block = 2 
  [../]
  [./creep2]
    type = LinearViscoelasticStressUpdate
     block = 2
  [../]
  [./logcreep]
    type = ConcreteLogarithmicCreepModel
    poissons_ratio = 0.2
    youngs_modulus = 12.00E9             
    recoverable_youngs_modulus = 12.00E9           
    recoverable_viscosity = 3.5E15            
    long_term_viscosity = 3.5E15            
    long_term_characteristic_time = 172800           
    temperature = temp
    activation_temperature = 5000
    humidity = h
    block = 2
  [../]
  [./damage2]
    type = IsotropicDamage
    cracking_release = exponential
    tensile_strength = 16E6                         
    fracture_energy=0.2E3             
    residual_fraction=0.0
    equivalent_strain_definition=modifiedvonmises
    use_old_damage= 'true'
    creep_damage_parameter=0.1
    block = 2
  [../]
  [./stress2]
    type = ComputeMultipleInelasticStress
    damage_model = damage2
    inelastic_models = 'creep2'
    block = 2
  [../]
[]

[UserObjects]
  [./soln1]
    type = SolutionUserObject
    mesh = PPT_ConcHum_vfpt3753_out.e
    system_variables = rh
  [../]
  [./soln2]
    type = SolutionUserObject
    mesh = PPT_ConcHum_vfpt3753_out.e
    system_variables = T
  [../]
  [./visco_update]
    type = LinearViscoelasticityManager
    viscoelastic_model = logcreep
    block = 2
  [../]
[]

[Contact]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart'
  petsc_options_value = 'asm lu 1 101'
  line_search = 'none'

  l_max_its = 50
  l_tol = 1.0e-5

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-5
  nl_max_its = 50

  dtmin = 86.4
  end_time = 13392000
  dt = 0.2
  [./TimeStepper]
    type = FunctionDT
    time_t =  ' 0 27648000'
    time_dt = '345600 345600'
  [../]
[]

[Outputs]
  file_base = ' Microstructure_with_radiation_out'
  exodus = true
  csv = false
[]
