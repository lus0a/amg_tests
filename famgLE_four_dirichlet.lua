--------------------------------------------------------------------------------
--[[!
-- \file apps/amg/famgLE_structured_unstretched.lua
--
-- test ob FAMG LE auf einheitsquadrat, mit eingespannt an einer Seite
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

bWriteStats = util.HasParamOption("-bWriteStats", "if true, write stats output")
bWriteMat = util.HasParamOption("-writeMatrices", "if true, write matrix output")

bAggressiveCoarsening = util.HasParamOption("-AC", "if true, use aggressive coarsening")
famgDelta = util.GetParamNumber("-delta", 0.2, "Interpolation quality F may not be worse than this (F < m_delta)")
famgTheta = util.GetParamNumber("-theta", 0.9, "with multiple parents, discard pairs with m_theta * F > min F.")

logtofile = util.GetParam("-logtofile", nil, "the file to write logs to")
bWCycle = util.HasParamOption("-Wcycle", "if true, use W-Cycle, otherwise: V-Cycle.")

iTestvectorSmooths = util.GetParamNumber("-iTestvectorSmooths", 0, "nr of testvector smooth for construction of P")
bUseJacobi = util.HasParamOption("-bUseJacobi", "if true, use jacobi as smoother. otherwise use gs")
iSmooths = util.GetParamNumber("-iSmooths", 2, "nr of pre/postsmooths")

gridName = util.GetParam("-grid", "amgtests/unit_square_quads_8x8_four_dirichlet_nodes.ugx")


--------------------------------------------------------------------------------------------------
util.CheckAndPrintHelp(
"Script for FAMG on Linear Elasticity\nMartin Rupp"
)
--------------------------------------------------------------------------------------------------

matrixWritePath=os.getenv("HOME").."/matrices/"

print("Writing matrices to "..matrixWritePath)
if bWriteMat and outdir == "" then
	outdir = matrixWritePath
end

dim = 2

-- choose algebra

InitUG(dim, AlgebraType("CPU", 2))

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
bLinearElasticity = true
approxSpace:add_fct("ux", "Lagrange", 1)          
approxSpace:add_fct("uy", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

elemDisc = LinearElastoPlasticity("ux, uy", "Inner")

-- set Elasticity-tensor "C" with material constants: 
-- first parameter: lambda ("1.lame constant")
-- second parameter: mu ("2.lame constant")
if util.HasParamOption("-bSmall") then
	elemDisc:set_hooke_elasticity_tensor(0.11074382,0.08019380)
else
	elemDisc:set_hooke_elasticity_tensor(11,8)
end

nu = util.GetParamNumber("-nu", 0.4)
elemDisc:set_hooke_elasticity_tensor_E_nu(20.6316, nu)

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

function DirichletValue2d(x, y, t)
	return 0.0
end

dirichletBND = DirichletBoundary()
function setBoundary(boundaryName)
dirichletBND:add("DirichletValue"..dim.."d", "ux", boundaryName)
dirichletBND:add("DirichletValue"..dim.."d", "uy", boundaryName)
end

setBoundary("Boundary")


-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)
-- domainDisc:add(neumannDisc)

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
jac = Jacobi(0.8)
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

base = LU()
base:set_minimum_for_sparse(0)
		
------- Create FAMG Solver ---------
amg = FAMGPreconditioner()
amg = util.SetupFAMGPreconditioner(base, presmoother, postsmoother)
amg:set_num_presmooth(iSmooths)
amg:set_num_postsmooth(iSmooths)

if bWCycle then
	amg:set_cycle_type(2)
else
	amg:set_cycle_type(1)
end

amg:set_write_f_values(true)

amg:set_aggressive_coarsening(true)

	
amg:set_max_nodes_for_base(maxBase)
amg:set_max_levels(maxLevels)


amg:set_testvector_smooths(iTestvectorSmooths)
amg:set_testvector_smoother(Jacobi(0.3))
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)


--amg:set_testvectors_LE_2d(1e3, 1e3)
amg:set_testvectors_linear_2d(1e3, 1e3, 1e3)
--amg:set_aggressive_coarsening(true)

amg:set_theta(famgTheta)

threePenalty = util.GetParamNumber("-threePenalty", 1.5)
maxParents = util.GetParamNumber("-maxParents", 2)
Hinv = util.GetParamBool("-Hinv", false)
delta = util.GetParamNumber("-delta", 0.1)

--amg:set_aggressive_coarsening(true)
amg:set_three_parents_penalty(threePenalty)
amg:set_max_parents(maxParents)
amg:set_H_as_inverse(Hinv)
amg:set_delta(delta)
--amg:set_delta(0.06, 0)

amg:set_one_init(true)
amg:set_injection_restriction( util.HasParamOption("-bUseInjection") )

if bWriteMat then
	amg:set_write_testvectors(true)
	amg:set_write_f_values(true)
	amg:set_matrix_write_path(matrixWritePath)
end	




vectorWriter = GridFunctionPositionProvider()
vectorWriter:set_reference_grid_function(u)
amg:set_position_provider(vectorWriter)

-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(50)
convCheck:set_minimum_defect(1e-15)
convCheck:set_reduction(1e-20)


if util.HasParamOption("-bYcycle") then
	convCheckY = ConvCheck()
	convCheckY:set_maximum_steps(4)
	convCheckY:set_minimum_defect(1e-30)
	convCheckY:set_reduction(0.1)
	convCheckY:set_verbose(false)
	
	amg:set_Y_check(convCheckY)
	amg:set_Y_cycle(40, 1e-30, 0.1)	
end




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

print(amg:config_string())
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)

amg:set_H_reduce_interpolation_nodes_parameter(0.5)	

-- 2. init solver for linear Operator

amg:set_one_init(true)

iSeed = 34

print("Apply solver: LINEAR SOLVER.")
	print("Init solver for operator.")
	linSolver = LinearSolver()
	if false then
		ilut = ILUTScalar(0.0)
		ilut:set_info(true)
		linSolver:set_preconditioner(ilut)
	else
		linSolver:set_preconditioner(amg)
	end
	linSolver:set_convergence_check(convCheck)
	linSolverA = AnalyzingSolver(linSolver)
	linSolverA:init(A)
	
	srand(iSeed)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	tBefore = os.clock()
	linSolver:apply_return_defect(u,b)	
	SaveVectorForConnectionViewer(u, "amg_u_linear.vec")
	tSolve = os.clock()-tBefore
	print("done")

	lastReduction = convCheck:defect()/convCheck:previous_defect()
	avgReduction = convCheck:avg_rate()
	steps = convCheck:step()
print("Apply solver: CG SOLVER.")
	cgSolver = CG()
	cgSolver:set_preconditioner(amg)
	cgSolver:set_convergence_check(convCheck)
	cgSolver:init(A)

	srand(iSeed)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	tBefore = os.clock()
	cgSolver:apply_return_defect(u,b)	
	tSolveCG = os.clock()-tBefore
	print("done")
	
	CGlastReduction = convCheck:defect()/convCheck:previous_defect()
	CGavgReduction = convCheck:avg_rate()
	CGsteps = convCheck:step()


if util.HasParamOption("-check") then
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	convCheck:set_maximum_steps(0)
	linSolver:apply_return_defect(u,b)
	--amg:check(u, b)
	amg:check()
end



-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 4 -stretchX 20  -delta 0.1 -Hinv 0 -maxParents 3  -threePenalty 1.4 --> 1.4

print(amg:level_information_str())

if GetProcessRank() == 0 then
	stats = {	
	{ "procs", GetNumProcesses() },
	--{ "numPreRefs", numPreRefs},
	{ "iGrid", iGrid },
	{ "numRefs", numRefs },
	{ "maxLevels", maxLevels},
	{ "stretchX", stretchX},
	{ "nu", nu},
	{ "delta", delta},
	{ "Hinv", Hinv },
	{ "maxParents", maxParents },
	{ "threePenalty", threePenalty },

	{ "steps", steps},
	{ "lastReduction", lastReduction},
	{ "avgReduction", avgReduction}, 
	
	{ "CGsteps", CGsteps},
	{ "CGlastReduction", CGlastReduction},
	{ "CGavgReduction", CGavgReduction}, 
	

	{ "c_A", amg:get_operator_complexity()},
	{ "c_G", amg:get_grid_complexity()},
	{ "used Levels", amg:get_used_levels()},
	{ "fineSize", amg:get_level_information(0):get_nr_of_nodes() },
	{ "coarseSize", amg:get_level_information(amg:get_used_levels()-1):get_nr_of_nodes() },
	 
	{ "date", os.date("y%Ym%md%d") },
	{ "SVN Revision", GetSVNRevision()},
	{"host",GetBuildHostname()},
	{"commandline", util.GetCommandLine() } } 
	
	util.printStats(stats)
	if bWriteStats  then	
		util.writeFileStats(stats, outdir.."stats.txt")
	end
end
--print(amg:config_string())
--print(elemDisc:config_string())	
print("\n")		
util.PrintArguments()
util.PrintIgnoredArguments()




--WriteProfileData("famgLE.pdxml")