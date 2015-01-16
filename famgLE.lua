--------------------------------------------------------------------------------
--[[!
-- \file apps/amg/famgLE.lua
--
-- \author Martin Rupp, Andreas Vogel
--
-- \brief Test Script for FAMG on Linear Elasticity
--
]]--
--------------------------------------------------------------------------------

-- ProfileLUA(true)
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util.lua")
ug_load_script("solver_util/setup_famg.lua")

-- command line parameters
--------------------------------------------------------------------------------------------------

paramOutputProcess = util.GetParamNumber("-outproc", 0, "Output processor")
numRefs = util.GetParamNumber("-numRefs", 1, "number of refinements")
numPreRefs = util.GetParamNumber("-numPreRefs", -1, "number of refinements before ")


maxBase = util.GetParamNumber("-maxBase", 2, "maximum number of nodes before base solver is used")
maxLevels = util.GetParamNumber("-maxLevels", 2, "maximum number of AMG levels")

RAepsilon = util.GetParamNumber("-RAepsilon", 1, "epsilon for rotated anisotropy")
RAalpha = util.GetParamNumber("-RAalpha", 0, "alpha (degrees) for rotated anisotropy")

epsx = util.GetParamNumber("-epsx", 1, "epsilon in X direction")
epsy = util.GetParamNumber("-epsy", 1, "epsilon in Y direction")


bWriteStats = util.HasParamOption("-bWriteStats", "if true, write stats output")
bWriteMat = util.HasParamOption("-writeMatrices", "if true, write matrix output")
bRSAMG = util.HasParamOption("-RSAMG", "if true, use RSAMG, otherwise FAMG") 

bAggressiveCoarsening = util.HasParamOption("-AC", "if true, use aggressive coarsening")
bExternalCoarsening = util.HasParamOption("-XC", "if true, use external coarsening in FAMG")
jacBaseSolver = util.HasParamOption("-jacBaseSolver", "if true, use Linear Jacobi as Base Solver")

famgDelta = util.GetParamNumber("-delta", 0.1, "Interpolation quality F may not be worse than this (F < m_delta)")
famgTheta = util.GetParamNumber("-theta", 0.9, "with multiple parents, discard pairs with m_theta * F > min F.")

minNodes = util.GetParamNumber("-minNodes", 10, "minimum number of nodes on one processor (for agglomeration)")
preferredNodes = util.GetParamNumber("-preferredNodes", 10, "preferred number of nodes on one processor (for agglomeration)")

outdir = util.GetParam("-outdir", "", "output directory")
b3d = util.HasParamOption("-3d", "if true, 3d, otherwise: 2d")
ninePoint = util.HasParamOption("-ninePoint", "for 2d: if true, use ninepoint, otherwise 5pt")
logtofile = util.GetParam("-logtofile", nil, "the file to write logs to")
bWCycle = util.HasParamOption("-Wcycle", "if true, use W-Cycle, otherwise: V-Cycle.")

iTestvectorSmooths = util.GetParamNumber("-iTestvectorSmooths", 2, "nr of testvector smooth for construction of P")
bUseJacobi = util.HasParamOption("-bUseJacobi", "if true, use jacobi as smoother. otherwise use gs")
iSmooths = util.GetParamNumber("-iSmooths", 2, "nr of pre/postsmooths")

numSys = util.GetParamNumber("-numSys", 1, "nr of systems")

-- gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
-- gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")
-- gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2.ugx")
gridName = util.GetParam("-grid", "unit_square_quads_8x8_four_dirichlet_nodes.ugx")
gridName = util.GetParam("-grid", "unit_square/neumann_stretched.ugx")
--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_16x16_stretched.ugx")
--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_16x16_scaled.ugx")
-- gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_dirichlet_bottom.ugx")
gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_16x16.ugx")
gridName = util.GetParam("-grid", "beam2d.ugx")
if util.HasParamOption("-sX") then
	gridName = util.GetParam("-grid", "beam2d_stretched_X.ugx")
elseif util.HasParamOption("-sY") then
	gridName = util.GetParam("-grid", "beam2d_stretched_Y.ugx")
end


--------------------------------------------------------------------------------------------------
util.CheckAndPrintHelp(
"Script for FAMG on Linear Elasticity\nMartin Rupp"
)
--------------------------------------------------------------------------------------------------


GetLogAssistant():set_output_process(paramOutputProcess)

matrixWritePath=os.getenv("HOME").."/matrices/"
print("Writing matrices to "..matrixWritePath)
if bWriteMat and outdir == "" then
	outdir = matrixWritePath
end

-- constants
if b3d == true then
	dim = 3
else
	dim = 2
end

iSeed = 23

nSystems = numSys

-- choose algebra
InitUG(dim, AlgebraType("CPU", 2))

--------------------------------
-- User Data Functions (begin)
--------------------------------
function ourDiffTensor2d(x, y, t)
	return	1, 0, 
			0, 1
end

function ourRhs2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
end

function ourDirichletBnd2d(x, y, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y)
end

function ourNeumannBnd(x, y, t)
	return 0
end

function ourDiffTensor3d(x, y, z, t)
	return	1, 0, 0,
			0, 1, 0,
			0, 0, 1
end


function ourRhs3d(x, y, z, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y) + math.sin(s*z))
end

function ourNeumannBnd3d(x, y, t)
	local s = 2*math.pi
	return -s*math.cos(s*x)
end

function ourDirichletBnd3d(x, y, z, t)
	local s = 2*math.pi
	return true, math.sin(s*x) + math.sin(s*y) + math.sin(s*z)
end
print ("Setting up1")
--------------------------------
-- User Data Functions (end)
--------------------------------

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets)

print ("Setting up2")
-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
bLinearElasticity = true
if bLinearElasticity then
	approxSpace:add_fct("ux", "Lagrange", 1)          
	approxSpace:add_fct("uy", "Lagrange", 1)
else  
	approxSpace:add_fct("c1", "Lagrange", 1)          
	approxSpace:add_fct("c2", "Lagrange", 1)
end
print ("Setting up3")
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")


--set_debug(debugID.RSAMG.Coarsening, 2)

-- Diffusion Tensor setup
--diffusionMatrix = LuaUserMatrix("ourDiffTensor"..dim.."d")
diffusionMatrix = ConstUserMatrix(1.0)

-- rhs setup
rhs = LuaUserNumber("ourRhs"..dim.."d")
--rhs = ConstUserNumber(0.0)


function ourNeumannBnd2d(x, y, t)
   return true, -1
end

-- dirichlet setup
dirichlet = LuaCondUserNumber("ourDirichletBnd"..dim.."d")
neumann = LuaCondUserNumber("ourNeumannBnd"..dim.."d")

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

if bLinearElasticity then	
		
	matLaw = HookeLaw()
	-- set Elasticity-tensor "C" with material constants: 
	-- first parameter: lambda ("1.lame constant")
	-- second parameter: mu ("2.lame constant")
	if util.HasParamOption("-bSmall") then
		matLaw:set_hooke_elasticity_tensor(0.11074382,0.08019380)
	else
		--matLaw:set_hooke_elasticity_tensor(110743.82,80193.80)
		matLaw:set_hooke_elasticity_tensor(11,8)
	end
	
	-- matLaw:set_hooke_elasticity_tensor(1.1,0.8)
	-- matLaw:set_hooke_elasticity_tensor(110000,80000)
	
	elemDisc = SmallStrainMechanics("ux, uy", "Inner")
	elemDisc:set_material_law(matLaw)
	-- set volume force "f":
	elemDisc:set_volume_forces(1) 
	
	-----------------------------------------------------------------
	--  Setup Dirichlet Boundary
	-----------------------------------------------------------------
	
	function DirichletValue2d(x, y, t)
		return 0.0
	end
	
	dirichletBND = DirichletBoundary()
	dirichletBND:add("DirichletValue"..dim.."d", "ux", "Boundary")
	dirichletBND:add("DirichletValue"..dim.."d", "uy", "Boundary")
	
	
	
	-- neumannZero = LuaCondUserNumber("ourNeumannBnd"..dim.."d")
	-- neumannDisc = NeumannBoundary("Inner")
	-- neumannDisc:add(neumann, "ux", "Boundary") -- not necessary when zero
	-- neumannDisc:add(neumann, "uy", "Boundary") -- not necessary when zero
	
	
	-------------------------------------------
	--  Setup Domain Discretization
	-------------------------------------------
	
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(elemDisc)
	domainDisc:add(dirichletBND)
	-- domainDisc:add(neumannDisc)
else
	
	diffusionMatrix = ConstUserMatrix(1.0)
	dirichlet = LuaCondUserNumber("ourDirichletBnd"..dim.."d")
	
	elemDisc = {}
	function Diffusion2d(x, y, t)
		return	1, 0, 
			0, 1
	end
	numSys = 2
	for i=1, numSys do
		elemDisc[i] = ConvectionDiffusion("c"..i, "Inner", "fv1")
		elemDisc[i]:set_diffusion("Diffusion2d")
	end
	dirichletBND = DirichletBoundary()
	for i=1, numSys do
		dirichletBND:add(dirichlet, "c"..i, "Boundary")
	end
	
	domainDisc = DomainDiscretization(approxSpace)
	for i=1, numSys do
		domainDisc:add(elemDisc[i])
	end
	domainDisc:add(dirichletBND)
end


-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- matrix and vectors
A = MatrixOperator()
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
b2 = GridFunction(approxSpace)

-- set initial value


-- init Operator and set dirichlet values in start iterate
print("Assemble A and b...")
tBefore = os.clock()
domainDisc:assemble_linear(A, b)
tAssemble = os.clock()-tBefore
print ("done")

-- write matrix for test purpose
--if bWriteMat then
	SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.vec")
--end

-- create algebraic Preconditioners

print ("create preconditioners... ")
jac = Jacobi()
jac:set_damp(0.8)

gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()


if bUseJacobi then
	presmoother = jac
	postsmoother = jac
else
	presmoother = gs
	postsmoother = gs
end




-- create Base Solver
baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(1000)
baseConvCheck:set_minimum_defect(1e-5)
baseConvCheck:set_reduction(1e-5)
baseConvCheck:set_verbose(false)

if false then
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)
else
	base = LU()
end

		
------- Create AMG Solver ---------
if true then
	amg = FAMGPreconditioner()
	amg = util.SetupFAMGPreconditioner(base, presmoother, postsmoother)
	amg:set_num_presmooth(2)
	amg:set_num_postsmooth(2)
	amg:set_write_f_values(true)
	
	if bExternalCoarsening then
		amg:set_external_coarsening(true)
		amg:set_parallel_coarsening(GetColorCoarsening())
		if dim == 2 then
			amg:set_strong_connection_external(0.2)
		else
			amg:set_strong_connection_external(0.1)
		end
	end
	amg:set_cycle_type(1)
	
	-- amg:set_delta(1e7)
	-- amg:set_delta(1e5)		
		
	amg:set_max_nodes_for_base(maxBase)
	amg:set_max_levels(maxLevels)
	-- amg:set_debug_level_calculate_parent_pairs(2)
	--amg:set_debug_level_precalculate_coarsening(4)
	-- amg:set_debug_level_external_coarsening(4)
	
	
	jac2 = Jacobi()
	jac2:set_damp(0.3)
	amg:set_testvector_smooths(1)
	amg:set_testvector_smoother(jac2)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
	
	
	
	if false then
		amg:set_testvector_from_matrix_rows()
	else
		-- amg:set_testvectors_LE(1e4)
		-- amg:set_testvectors_LE(1e7) -- gg
		-- amg:set_testvectors_LE_2d(1e3, 1e3)
		amg:set_testvectors_LE_2d(1e3, 1e3)
		bWriteMat = true		
	end	
	
	if false then
		-- neumann.ugx
		amg:set_delta(0.05)
		amg:set_theta(0.9)
		amg:set_testvectors_LE(1e3)
	else
		-- neumann_stretched
		--amg:set_delta(0.26)
		amg:set_delta(0.4)
		amg:set_theta(0.9)
		amg:set_testvectors_LE_2d(1e3, 1e3)
		--amg:set_testvector_from_matrix_rows()
	end
	--amg:set_delta(1.4)
	
	--amg:set_aggressive_coarsening(true)
	case = 2
	if case == 1 then
		amg:set_H_as_inverse(true)
		amg:set_delta(2)
	elseif case == 2 then
		amg:set_H_as_inverse(false)
		amg:set_aggressive_coarsening(false)
		amg:set_delta(0.3)
	else
		amg:set_H_as_inverse(true)
		amg:set_aggressive_coarsening(true)
		--amg:set_delta(1.4)
		amg:set_delta(1.4)
	
	end
	
	
else

	print ("create AMG... ")
	amg = RSAMGPreconditioner()
	-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
	amg:set_parallel_coarsening(GetColorCoarsening())
	-- amg:set_parallel_coarsening(GetRS3Coarsening())
	-- amg:set_parallel_coarsening(GetSimpleParallelCoarsening())
	if bAggressiveCoarsening then
		amg:enable_aggressive_coarsening_A(2)
	end
	amg:set_presmoother(gs)
	amg:set_postsmoother(gs)
	amg:set_base_solver(base)
end
amg:set_one_init(true)

if bWriteMat then
	vectorWriter = GridFunctionPositionProvider()
	vectorWriter:set_reference_grid_function(u)
	amg:set_position_provider(vectorWriter)

	amg:set_write_testvectors(true)
	amg:set_write_f_values(true)
	amg:set_matrix_write_path(matrixWritePath)
end	

		

amg:set_injection_restriction(false)
-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(60)
convCheck:set_minimum_defect(1e-15)
convCheck:set_reduction(1e-20)

print("done.")
-- create Linear Solver

-- create GMG solver --




-- Apply Solver

b2 = GridFunction(approxSpace)
b2:assign(b)

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator


if true then
	dbgWriter = GridFunctionDebugWriter(approxSpace)
	dbgWriter:set_base_dir(matrixWritePath)
	dbgWriter:set_vtk_output(false)
	transfer = StdTransfer()
	transfer:set_debug(dbgWriter)
	
	
	-- Gemoetric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_transfer(transfer)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0) -- 1
	gmg:set_base_solver(base)
	gmg:set_smoother(gs)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)
	gmgSolver = LinearSolver()
	gmgSolver:set_preconditioner(gmg)
	gmgCGSolver = CG()
	gmgCGSolver:set_preconditioner(gmg)
	
	-- create Convergence Check
	gmgConvCheck = ConvCheck()
	gmgConvCheck:set_maximum_steps(50)
	gmgConvCheck:set_minimum_defect(1e-24)
	gmgConvCheck:set_reduction(1e-24)
	
	gmgSolver:set_convergence_check(gmgConvCheck)
	gmgCGSolver:set_convergence_check(gmgConvCheck)
	
	
	
	u = GridFunction(approxSpace)
	b = GridFunction(approxSpace)
	u:set_random(-1.0, 1.0)
	
	
	
	linOp = AssembledLinearOperator()
	linOp:set_discretization(domainDisc)

	-- 1. init operator
	print("GMG/CG");
	if AssembleLinearOperatorRhsAndSolution(linOp, u, b) == false then 
		print("Could not assemble operator"); exit(); 
	end
	
	b2 = GridFunction(approxSpace)
	b2:assign(b)
	
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	domainDisc:adjust_solution(u)
	
		-- 2. apply solver
	if false then
		print("Apply solver.")
		if ApplyLinearSolver(linOp, u, b, gmgCGSolver) == false then
			print("Could not apply linear solver.");
		end
		
		print("done");
	end
	uEx = GridFunction(approxSpace)
	uEx:assign(u)
	
		-- 1. init operator
	print("GMG");
	b:assign(b2)
	
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	domainDisc:adjust_solution(u)
		-- 2. apply solver
	print("Apply solver.")
	--gmgConvCheck:set_maximum_steps()
	
	
	--gmg:set_debug(dbgWriter)
	
	
	if ApplyLinearSolver(linOp, u, b, gmgSolver) == false then
		print("Could not apply linear solver.");
	end
	SaveVectorDiffForConnectionViewer(u, uEx, "gmg_diff.vec")
	print("done");
	
	
	SaveVectorForConnectionViewer(u, "gmg_u.vec")
	--SaveVectorForConnectionViewer(u, A, "gmg_u.vec")	
	SaveVectorForConnectionViewer(b, A, "gmg_b.vec")
end
	print(amg:tostring())
	amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
	amg:set_max_parents(3)
	amg:set_H_reduce_interpolation_nodes_parameter(0.5)	

-- 2. init solver for linear Operator
print("Init solver for operator.")
	linSolver = LinearSolver()
	linSolver:set_preconditioner(amg)
	linSolver:set_convergence_check(convCheck)
	linSolver:init(A)

	

	print("Apply solver: LINEAR SOLVER.")
	srand(iSeed)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	tBefore = os.clock()
	linSolver:apply_return_defect(u,b)	
	SaveVectorForConnectionViewer(u, "amg_u20.vec")
	tSolve = os.clock()-tBefore
	print("done")
	
	
	lastReduction = convCheck:defect()/convCheck:previous_defect()
	steps = convCheck:step()
	solution = GridFunction(approxSpace)
	solution:assign(u)
	SaveVectorDiffForConnectionViewer(u, uEx, "amg_diff.vec")
	SaveVectorForConnectionViewer(u, A, "amg_u.vec")
	SaveVectorForConnectionViewer(u, "amg_u2.vec")		
	SaveVectorForConnectionViewer(b, A, "amg_b.vec")
	amg:check_interpolation(u, "u")
	
exit()
---------
if true then
	linSolver = CG()
	linSolver:set_preconditioner(amg)
	linSolver:set_convergence_check(convCheck)
	linSolver:init(A)
	
	print("Apply solver: CG.")
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	SaveVectorForConnectionViewer(u, A, "rand_u.vec")
	domainDisc:adjust_solution(u)
	b:assign(b2)
	tBefore = os.clock()
	linSolver:apply_return_defect(u,b)
	tSolveCG = os.clock()-tBefore
	print("done")
	
	lastReductionCG = convCheck:defect()/convCheck:previous_defect()
	stepsCG = convCheck:step()
	
	SaveVectorDiffForConnectionViewer(u, uEx, "amg_cg_diff.vec")
	SaveVectorForConnectionViewer(u, A, "amg_cg_u.vec")	
	SaveVectorForConnectionViewer(b, A, "amg_cg_b.vec")
end
---------
if util.HasParamOption("-bCheck") then
	print("CHECKS:")
	
	linSolver = LinearSolver()	
	linSolver:set_preconditioner(amg)
	linSolver:set_convergence_check(convCheck)
	linSolver:init(A)
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	domainDisc:adjust_solution(u)
	b:assign(b2)

	convCheck:set_maximum_steps(2)
	linSolver:apply_return_defect(u,b)
	
	if bWriteMat then
		SaveVectorForConnectionViewer(b, A, "b.vec")
		SaveVectorForConnectionViewer(solution, A, "solution.vec")
		SaveVectorForConnectionViewer(u, solution, A, "u-solution.vec")
		amg:write_interfaces()
	end
	print("amg:check(u, b)")
	amg:check(u,b)	
	--GetLogAssistant():set_output_process(util.GetParamNumber("-outproc", 0))
	
	--print("amg:check_testvector()")
	-- amg:check_testvector()
	--print("amg:check_fsmoothing()")
	-- amg:check_fsmoothing()
end	
	
---------



if ProcRank() == 0 then
	stats = {
	
	{ "procs", NumProcs() },
	--{ "numPreRefs", numPreRefs},
	{ "numRefs", numRefs },
	{ "RAepsilon", RAepsilon},
	{ "dim", dim},
	{ "gridName", gridName},
	{ "ndofs", amg:get_level_information(0):get_nr_of_nodes() },
	{ "tSetupAmg [s]", amg:get_timing_whole_setup_ms()/1000},
	{ "main [ms]", GetProfileNode("main"):get_avg_total_time_ms()},
	{ "perform_refinement [ms]", GetProfileNode("perform_refinement"):get_avg_total_time_ms()},
	{ "c_create_AMG_level [ms]", GetProfileNode("c_create_AMG_level"):get_avg_total_time_ms()},
	{ "init_fsmoothing [ms]", GetProfileNode("init_fsmoothing"):get_avg_total_time_ms()},
	{ "amg_createDirectSolver [ms]", GetProfileNode("amg_createDirectSolver"):get_avg_total_time_ms()},
	
	
	{ "XC", bExternalCoarsening},	
	{ "AC", bAggressiveCoarsening},
	{ "c_A", amg:get_operator_complexity()},
	{ "c_G", amg:get_grid_complexity()},
	{ "used Levels", amg:get_used_levels()},
	 
	{ "tSolve [s]", tSolve},
	{ "steps", steps},
	{ "lastReduction", lastReduction},
	
	{ "tSolveCG [s]", tSolveCG},
	{ "stepsCG", stepsCG},
	{ "lastReductionCG", lastReductionCG},
	
	{ "tGrid [s]",  tGrid},	
	{ "tAssemble [s]", tAssemble},
	
	{ "maxBase", maxBase},
	{ "date", os.date("y%Ym%md%d") },
	{ "SVN Revision", GetSVNRevision()},
	{"host",GetBuildHostname()},
	{"filename", filename},
	{"commandline", util.GetCommandLine() } } 
	
	util.printStats(stats)
	if bWriteStats  then	
		util.writeFileStats(stats, outdir.."stats.txt")
	end

end
	
if not bCheck and bWriteStats then
		-- ..os.date("y%Ym%md%d").. -- table.insert
	
	if GetProfilerAvailable() == true then
		create_levelPN = GetProfileNode("c_create_AMG_level")
		
		if true then
			print(create_levelPN:call_tree())
			print(create_levelPN:child_self_time_sorted())
		end
		if true then			
			to100 = create_levelPN:get_avg_total_time_ms()
			PrintParallelProfileNode("create_OL2_matrix")
			PrintParallelProfileNode("CalculateTestvector")
			PrintParallelProfileNode("CreateSymmConnectivityGraph")
			PrintParallelProfileNode("calculate_all_possible_parent_pairs")
			PrintParallelProfileNode("color_process_graph")
			PrintParallelProfileNode("FAMG_recv_coarsening_communicate")
			PrintParallelProfileNode("update_rating")
			PrintParallelProfileNode("precalculate_coarsening")
			PrintParallelProfileNode("send_coarsening_data_to_processes_with_higher_color")
			PrintParallelProfileNode("communicate_prolongation")
			PrintParallelProfileNode("create_new_index")
			PrintParallelProfileNode("create_parent_index")
			PrintParallelProfileNode("create_interfaces")
			PrintParallelProfileNode("create_fine_marks")
			PrintParallelProfileNode("create_galerkin_product")
			PrintParallelProfileNode("CalculateNextTestvector")
		end
		
	else
		print("Profiler not available.")
	end	
end

print("\n")		
util.PrintArguments()
util.PrintIgnoredArguments()
WriteProfileData("famgLE.pdxml")