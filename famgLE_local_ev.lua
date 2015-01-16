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

InitSignals()
	
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

gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_4bnd.ugx")
--gridName = util.GetParam("-grid", "amgtests/unstructured_33k.ugx")
--gridName = util.GetParam("-grid", "amgtests/unstructured_3k.ugx")


--------------------------------------------------------------------------------------------------
util.CheckAndPrintHelp(
"Script for FAMG on Linear Elasticity\nMartin Rupp"
)
--------------------------------------------------------------------------------------------------

if bWriteStats then
	outdir=os.getenv("HOME").."/results/famgLE/"
	commandLineFilename = util.GetUniqueFilenameFromCommandLine()
	thePath = outdir.."famgLE"..commandLineFilename.."/"
	CreateDirectory(thePath)
	
	GetLogAssistant():enable_file_output(false, "")
	GetLogAssistant():enable_file_output(true, thePath.."log.txt")
else
	thePath=os.getenv("HOME").."/matrices/"
end
matrixWritePath=thePath


print("Writing matrices to "..matrixWritePath)
if bWriteMat and outdir == "" then
	outdir = matrixWritePath
end

dim = 2

-- choose algebra

InitUG(dim, AlgebraType("CPU", 2))

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", "Bottom"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, 0, neededSubsets)

stretchX = util.GetParamNumber("-stretchX", 1.0)
stretchY = util.GetParamNumber("-stretchY", 1.0) 
ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))
--d = math.pow(0.3, numRefs)
d = math.pow(0.25, numRefs)
d = math.pow(0.1, numRefs)
d = math.pow(0.25, numRefs)
print("randomize ="..d)
--RandomizeDomain(dom, d, d, d)
print("stretchX = "..stretchX..", stretchY = "..stretchY)

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
	--elemDisc:set_hooke_elasticity_tensor(0.11074382,0.08019380)
else
	--elemDisc:set_hooke_elasticity_tensor(11,8)
end

nu = util.GetParamNumber("-nu", 0.4)
--elemDisc:set_hooke_elasticity_tensor_E_nu(20.6316, nu)
elemDisc:set_hooke_elasticity_tensor_E_nu(0.001*20.6316, nu)

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

Boundaries = {"Bottom", "Top", "Left", "Right"}
bnd = util.GetParamNumber("-bnd", 1)
setBoundary("Bottom")
if bnd > 1 then
setBoundary("Top")
if bnd > 2 then
setBoundary("Left")
if bnd > 3 then
setBoundary("Right")
end end end

OrderCuthillMcKee(approxSpace, false);

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
bWriteMat = true
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

--postsmoother = BlockGaussSeidel(1)
postsmoother = GaussSeidel()

presmoother = postsmoother




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



	
amg:set_max_nodes_for_base(maxBase)
amg:set_max_levels(maxLevels)


amg:set_testvector_smooths(iTestvectorSmooths)
amg:set_testvector_smoother(Jacobi(0.3))
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)

if util.HasParamOption("-XC") then
amg:set_external_coarsening(true)
end



if util.HasParamOption("-tvLin") then
tvDesc = "tvLin"
amg:set_testvectors_linear_2d(1e3, 1e3, 1e3)
elseif util.HasParamOption("-tvL") then
tvDesc = "tvLComp"
amg:set_testvectors_linear_component()
elseif util.HasParamOption("-tvLE") then
tvDesc = "tvLE"
amg:set_testvectors_LE_2d(1e3, 1e3)
else
tvDesc = "tvConst"
end
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
convCheck:set_maximum_steps(4)
convCheck:set_minimum_defect(1e-15)
convCheck:set_reduction(1e-20)
convCheck:set_verbose(true)	


if util.HasParamOption("-bYcycle") then
	convCheckY = ConvCheck()
	convCheckY:set_maximum_steps(4)
	convCheckY:set_minimum_defect(1e-30)
	convCheckY:set_reduction(0.1)
	convCheckY:set_verbose(false)
	
	--amg:set_Y_check(convCheckY)
	amg:set_Y_cycle(40, 1e-30, 0.2)	
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




--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")


-- choose precond
print("Init solver for operator.")

results = {}
		
function bla(i)


linSolver = LinearSolver()

sm = SymmetricBlockGaussSeidelIterative(i, 1)

--sm = SymmetricGaussSeidel()
	linSolver:set_preconditioner(sm)	
	linSolver:set_convergence_check(convCheck)
	linSolver:init(A)		

convCheck:set_maximum_steps(20)
	b:assign(b2)
	SaveVectorForConnectionViewer(b, "bb"..i..".vec")
	--srand()
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	SaveVectorForConnectionViewer(u, "u0"..i..".vec")
	tBefore = os.clock()	
	linSolver:apply_return_defect(u,b)
	tSolve = os.clock()-tBefore
	print("took "..tSolve.." s")
	SaveVectorForConnectionViewer(u, "smooth"..i..".vec")
	SaveVectorForConnectionViewer(b, "b"..i..".vec")
end	
	
	
for i=1,3 do
bla(i)
end