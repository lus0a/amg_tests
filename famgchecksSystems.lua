--------------------------------------------------------------------------------
--[[!
-- \file apps/amg/famgchecksSystems.lua
--
-- \author Martin Rupp, Andreas Vogel
--
-- \brief Test Script for FAMG for Systems
--
]]--
--------------------------------------------------------------------------------

-- ProfileLUA(true)
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util.lua")

-- command line parameters
--------------------------------------------------------------------------------------------------

paramOutputProcess = util.GetParamNumber("-outproc", 0, "Output processor")
numRefs = util.GetParamNumber("-numRefs", 1, "number of refinements")
numPreRefs = util.GetParamNumber("-numPreRefs", -1, "number of refinements before ")


maxBase = util.GetParamNumber("-maxBase", 600, "maximum number of nodes before base solver is used")
maxLevels = util.GetParamNumber("-maxLevels", 30, "maximum number of AMG levels")

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

famgDelta = util.GetParamNumber("-delta", 0.1, "\"Interpolation quality\" F may not be worse than this (F < m_delta)")
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

gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")


--------------------------------------------------------------------------------------------------
util.CheckAndPrintHelp(
"FAMG Script by Martin Rupp"
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
InitUG(dim, AlgebraType("CPU", 1))

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
dom = util.CreateAndDistributeDomain(gridName, numRefs-1, numPreRefs, neededSubsets)

print ("Setting up2")
-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
for i=1, nSystems do
	approxSpace:add_fct("c"..i, "Lagrange", 1)
end
print ("Setting up3")
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")


-- Diffusion Tensor setup
--diffusionMatrix = LuaUserMatrix("ourDiffTensor"..dim.."d")
diffusionMatrix = ConstUserMatrix(1.0)

-- rhs setup
rhs = LuaUserNumber("ourRhs"..dim.."d")
--rhs = ConstUserNumber(0.0)

-- dirichlet setup
dirichlet = LuaCondUserNumber("ourDirichletBnd"..dim.."d")

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------
elemDisc = {}
upwind = {}
for i=1, nSystems do
	upwind[i] = FullUpwind2d()
	elemDisc[i] = ConvectionDiffusion("c"..i, "Inner")
	elemDisc[i]:set_disc_scheme("fv1")
	elemDisc[i]:set_upwind(upwind[i])
	elemDisc[i]:set_diffusion(diffusionMatrix)
	elemDisc[i]:set_source(rhs)
end

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

dirichletBND = DirichletBoundary()
for i=1, nSystems do
	dirichletBND:add(dirichlet, "c"..i, "Boundary")
end

-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
for i=1, nSystems do
	domainDisc:add(elemDisc[i])
end
domainDisc:add(dirichletBND)


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
if bWriteMat then
	SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.vec")
end

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
	presmoother = sgs
	postsmoother = sgs
end




-- create Base Solver
baseConvCheck = StandardConvergenceCheck()
baseConvCheck:set_maximum_steps(1000)
baseConvCheck:set_minimum_defect(1e-16)
baseConvCheck:set_reduction(1e-16)
baseConvCheck:set_verbose(false)

if false then
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)
else
	base = LU()
end


-- create AMG ---
-----------------

if bRSAMG == false then
	print ("create FAMG... ")
	-- Testvectors for FAMG ---
	--------------------------	
	function CreateAMGTestvector(gridfunction, luaCallbackName, dim)
		local amgTestvector;
		amgTestvector = GridFunctionVectorWriter()
		amgTestvector:set_reference_grid_function(gridfunction)
		amgTestvector:set_user_data(LuaUserNumber(luaCallbackName))
		return amgTestvector	
	end
		
	function CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
		local amgDirichlet0 = GridFunctionVectorWriterDirichlet0()
		amgDirichlet0:init(dirichletBND, approxSpace)
		return amgDirichlet0
	end
	
	amg = FAMGPreconditioner()
	
	amg:set_delta(famgDelta)
	amg:set_theta(famgTheta)
	amg:set_aggressive_coarsening(bAggressiveCoarsening)
		
	-- add testvector which is 1 everywhere and only 0 on the dirichlet Boundary.
	testvectorwriter = CreateAMGTestvectorDirichlet0(dirichletBND, approxSpace)
	testvector = GridFunction(approxSpace)
	testvectorwriter:update(testvector)	
	amg:add_testvector(testvectorwriter, 1.0)


	jac2 = Jacobi()
	jac2:set_damp(0.8)
	amg:set_testvector_smooths(iTestvectorSmooths)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
	amg:set_testvector_smoother(jac2)
		
	if bWriteMat then
		amg:set_write_testvectors(true)
		amg:set_write_f_values(true)
	end
	
	if bExternalCoarsening then
		amg:set_external_coarsening(true)
		amg:set_parallel_coarsening(GetColorCoarsening())
		if dim == 2 then
			amg:set_strong_connection_external(0.2)
		else
			amg:set_strong_connection_external(0.1)
		end
		-- amg:set_parallel_coarsening(GetFullSubdomainBlockingCoarsening())
		-- amg:set_parallel_coarsening(GetRS3Coarsening())
	end
	
	
	-- amg:set_debug_level_overlap(4, 4)
	-- amg:set_use_precalculate(false)
	-- amg:set_debug_level_get_ratings(4)
	-- amg:set_debug_level_coloring(4)
	-- amg:set_debug_level_communicate_prolongation(4)
	-- amg:set_debug_level_overlap(4,4)
	-- amg:set_debug_level_calculate_parent_pairs(6)
	-- amg:set_debug_level_calculate_parent_pairs(6)
	-- amg:set_debug_level_precalculate_coarsening(4)
	-- amg:set_galerkin_truncation(1e-6)
	-- amg:set_H_reduce_interpolation_nodes_parameter(0.1)
	amg:set_galerkin_truncation(1e-9)
	amg:set_H_reduce_interpolation_nodes_parameter(0.0) --1)
	
	if bExternalCoarsening then
		amg:set_prereduce_A_parameter(0.1)
	else
		amg:set_prereduce_A_parameter(0.0)	
	end
	
	
	amg:set_nr_of_preiterations_at_check(10)
	amg:set_preiterations_mimum_defect_at_check(1e-10)	
	-- amg:set_write_coarsening()
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
end


vectorWriter = GridFunctionPositionProvider()
vectorWriter:set_reference_grid_function(u)
amg:set_position_provider(vectorWriter)
if bWriteMat then
amg:set_matrix_write_path(matrixWritePath)
end



amg:set_num_presmooth(iSmooths)
amg:set_num_postsmooth(iSmooths)
if bWCycle then
	amg:set_cycle_type(2)
else
	amg:set_cycle_type(1)
end
-- amg:set_Y_cycle(2, 1e-1, 1e-12)
amg:set_presmoother(presmoother)
amg:set_postsmoother(postsmoother)
amg:set_base_solver(base)
amg:set_max_levels(maxLevels)

amg:set_min_nodes_on_one_processor(minNodes)
amg:set_preferred_nodes_on_one_processor(preferredNodes)
-- amg:set_preferred_nodes_on_one_processor(100000)

amg:set_max_nodes_for_base(maxBase)
amg:set_max_fill_before_base(0.5)
amg:set_fsmoothing(true)
amg:set_prolongation_truncation(0)
amg:tostring()


-- create Convergence Check
convCheck = StandardConvergenceCheck()
convCheck:set_maximum_steps(20)
convCheck:set_minimum_defect(1e-30)
convCheck:set_reduction(1e-30)

print("done.")
-- create Linear Solver


-- Apply Solver

b2 = GridFunction(approxSpace)
b2:assign(b)

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator


amg:set_one_init(true)

-- 2. init solver for linear Operator
print("Init solver for operator.")
	linSolver = LinearSolver()
	linSolver:set_preconditioner(amg)
	linSolver:set_convergence_check(convCheck)
	linSolver:init(A)

	print("Apply solver: LINEAR SOLVER.")
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	domainDisc:adjust_solution(u)
	tBefore = os.clock()
	linSolver:apply_return_defect(u,b)	
	tSolve = os.clock()-tBefore
	print("done")
	
	lastReduction = convCheck:defect()/convCheck:previous_defect()
	steps = convCheck:step()
	solution = GridFunction(approxSpace)
	solution:assign(u)

---------
if false then
	linSolver = CG()
	linSolver:set_preconditioner(amg)
	linSolver:set_convergence_check(convCheck)
	linSolver:init(A)
	
	print("Apply solver: CG.")
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	domainDisc:adjust_solution(u)
	b:assign(b2)
	tBefore = os.clock()
	linSolver:apply_return_defect(u,b)
	tSolveCG = os.clock()-tBefore
	print("done")
	
	lastReductionCG = convCheck:defect()/convCheck:previous_defect()
	stepsCG = convCheck:step()
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
	--amg:check(u,b)	
	--GetLogAssistant():set_output_process(util.GetParamNumber("-outproc", 0))
	
	--print("amg:check_testvector()")
	amg:check_testvector()
	--print("amg:check_fsmoothing()")
	-- amg:check_fsmoothing()
end	
	
---------



if GetProcessRank() == 0 then
	stats = {
	
	{ "procs", GetNumProcesses() },
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