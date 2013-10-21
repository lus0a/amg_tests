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

bWriteStats = util.HasParamOption("-bWriteStats", "if true, write stats output")
bWriteMat = util.HasParamOption("-writeMatrices", "if true, write matrix output")

bAggressiveCoarsening = util.HasParamOption("-AC", "if true, use aggressive coarsening")
famgDelta = util.GetParamNumber("-delta", 0.1, "Interpolation quality F may not be worse than this (F < m_delta)")
famgTheta = util.GetParamNumber("-theta", 0.9, "with multiple parents, discard pairs with m_theta * F > min F.")

outdir = util.GetParam("-outdir", "", "output directory")
logtofile = util.GetParam("-logtofile", nil, "the file to write logs to")
bWCycle = util.HasParamOption("-Wcycle", "if true, use W-Cycle, otherwise: V-Cycle.")

iTestvectorSmooths = util.GetParamNumber("-iTestvectorSmooths", 2, "nr of testvector smooth for construction of P")
bUseJacobi = util.HasParamOption("-bUseJacobi", "if true, use jacobi as smoother. otherwise use gs")
iSmooths = util.GetParamNumber("-iSmooths", 2, "nr of pre/postsmooths")

gridName = util.GetParam("-grid", "beam2d.ugx")

stretchX = util.GetParamNumber("-stretchX", 1.0)
stretchY = util.GetParamNumber("-stretchY", 1.0) 


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

dim = 2
iSeed = 23

-- choose algebra

InitUG(dim, AlgebraType("CPU", 2))

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets)

-- scale domain
ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))

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

elemDisc:set_volume_forces(0) 

-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

function DirichletValue2d(x, y, t)
	return 0.0
end

dirichletBND = DirichletBoundary()
dirichletBND:add("DirichletValue"..dim.."d", "ux", "Boundary")
dirichletBND:add("DirichletValue"..dim.."d", "uy", "Boundary")

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
-- amg:set_debug_level_precalculate_coarsening(4)
-- amg:set_debug_level_external_coarsening(4)


jac2 = Jacobi()
jac2:set_damp(0.3)
amg:set_testvector_smooths(0)
amg:set_testvector_smoother(jac2)
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)


 -- always calc 3
--amg:set_max_parents(4) -- calc 3 if 2 are < 0.5 delta



-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 0.1  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.5 ---> 0.032
-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 0.1  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.5 ---> 0.044

-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 0.5  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.5 ---> 0.21
-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 0.5  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.5 ---> 0.21

-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.5 ---> 0.24
-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.0 ---> 0.22
-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 2.0 ---> 0.17
-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.3 -Hinv 0 -maxParents 2 -threePenalty 2.0 ---> 0.6
-- ugshell -ex amg/famgLEcheck.lua -numRefs 3 -maxLevels 3 -stretchX 20  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 ---> 0.029

amg:set_testvectors_LE_2d(1e3, 1e3)

amg:set_theta(famgTheta)

threePenalty = util.GetParamNumber("-threePenalty", 1.5)
maxParents = util.GetParamNumber("-maxParents", 2)
Hinv = util.GetParamBool("-Hinv", false)
delta = util.GetParamNumber("-delta", 0.1)

--amg:set_aggressive_coarsening(true)
amg:set_three_parents_penalty(threePenalty)
	amg:set_max_parents(maxParents)
	amg:set_H_as_inverse(Hinv)
	--amg:set_aggressive_coarsening(true)
	amg:set_delta(delta)

--[[
case = 2
if case == 1 then
	amg:set_H_as_inverse(true)
	amg:set_delta(2)
elseif case == 2 then
	amg:set_H_as_inverse(false)
	amg:set_aggressive_coarsening(true)
	amg:set_delta(0.2)
else
	amg:set_H_as_inverse(true)
	amg:set_aggressive_coarsening(true)
	--amg:set_delta(1.4)
	amg:set_delta(1.4)

end
]]--

amg:set_one_init(true)


	vectorWriter = GridFunctionPositionProvider()
	vectorWriter:set_reference_grid_function(u)
	amg:set_position_provider(vectorWriter)

--if bWriteMat then
	amg:set_write_testvectors(true)
	amg:set_write_f_values(true)
	amg:set_matrix_write_path(matrixWritePath)
--end	

		

amg:set_injection_restriction(false)
-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(30)
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

print(amg:tostring())
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)

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


if util.HasParamOption("-check") then
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	convCheck:set_maximum_steps(4)
	linSolver:apply_return_defect(u,b)
	amg:check(u, b)
end

lastReduction = convCheck:defect()/convCheck:previous_defect()
steps = convCheck:step()

-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 4 -stretchX 20  -delta 0.1 -Hinv 0 -maxParents 3  -threePenalty 1.4 --> 1.4

if GetProcessRank() == 0 then
	stats = {	
	{ "procs", GetNumProcesses() },
	--{ "numPreRefs", numPreRefs},
	{ "numRefs", numRefs },
	{ "maxLevels", maxLevels},
	{ "stretchX", stretchX},
	{ "delta", delta},
	{ "Hinv", Hinv },
	{ "maxParents", maxParents },

	{ "steps", steps},
	{ "lastReduction", lastReduction},

	{ "c_A", amg:get_operator_complexity()},
	{ "c_G", amg:get_grid_complexity()},
	{ "used Levels", amg:get_used_levels()},
	 
	{ "date", os.date("y%Ym%md%d") },
	{ "SVN Revision", GetSVNRevision()},
	{"host",GetBuildHostname()},
	{"commandline", util.GetCommandLine() } } 
	
	util.printStats(stats)
	if bWriteStats  then	
		util.writeFileStats(stats, outdir.."stats.txt")
	end
end
	
print("\n")		
util.PrintArguments()
util.PrintIgnoredArguments()
--WriteProfileData("famgLE.pdxml")