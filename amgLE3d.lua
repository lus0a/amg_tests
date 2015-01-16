--------------------------------------------------------------------------------
--[[!
-- \file cembalo/cembaloSolve.lua
-- \author Martin Rupp
-- \brief Solver script for the Cembalo
--
]]--
--------------------------------------------------------------------------------

PrintBuildConfiguration()

ug_load_script("ug_util.lua")
ug_load_script("solver_util/setup_rsamg.lua")
ug_load_script("solver_util/setup_famg.lua")

--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------

-- command line parameters
--------------------------------------------------------------------------------------------------
matrixWritePath=os.getenv("HOME").."/matrices/"

paramOutputProcess = util.GetParamNumber("-outproc", 0, "Output processor")
numRefs = util.GetParamNumber("-numRefs", 1, "number of refinements")
numPreRefs = util.GetParamNumber("-numPreRefs", -1, "number of refinements before ")


maxBase = util.GetParamNumber("-maxBase", 10, "maximum number of nodes before base solver is used")
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

-- Several definitions which can be changed by command line parameters

dim = 3 -- dim is always 3

-- CEMBALO PARAMETERS --



-- grid
gridName = util.GetParam("-grid", "cembalo/cembalo_3d_test2.ugx")

-- rho = density
rho = util.GetParamNumber("-rho", 0.00050, "density")

-- orthrope linear elasticity tensor
C11 = 7341919.0868408
C12 = 4908133.7541864
C13 = 2729895.8814608
C22 = 163359543.45898
C23 = 4287702.8720363
C33 = 7910724.2195105
C44 = 6167704.4025157
C55 = 363209.25925926
C66 = 7602054.2635659

if false then
	C11 =   7.3
	C12 =   4.9
	C13 =   2.7
	C22 = 163.4
	C23 =   4.3
	C33 =   7.9
	C44 =   6.2
	C55 =   0.4
	C66 =   7.6
end

-- eigensolver parameters:

-- refinements:
numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
numRefs    = util.GetParamNumber("-numRefs",    0, "number of refinements")


if numPreRefs > numRefs then
	print("It must be choosen: numPreRefs <= numRefs");
	exit();
end

--! way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "bisect", "way the domain will be distributed: [grid2d | bisect | metis]")

--! amount of output
--! set to 0 i.e. for time measurements,
--! >= 1 for writing matrix files etc.
verbosity = util.GetParamNumber("-verb", 0, "verbosity. 0 for time measurements")		

activateDbgWriter = 0	  
--! set to 0 i.e. for time measurements,
--! >= 1 for debug output: call 'set_debug(dbgWriter)'
--! for the main solver ('gmg')
activateDbgWriter = util.GetParamNumber("-dbgw", 0, "debug writer: 0 for no output.")


util.CheckAndPrintHelp("Linear Elasticity Eigenvalue\nMartin Rupp, Rossitza Piperkova");

-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)
print("  C11 = "..C11)
print("  C12 = "..C12)
print("  C13 = "..C13)
print("  C22 = "..C22)
print("  C23 = "..C23)
print("  C33 = "..C33)
print("  C44 = "..C44)
print("  C55 = "..C55)
print("  C66 = "..C66)
print("  rho = "..rho)

print("    verb (verbosity)         = " .. verbosity)
print("    dbgw (activateDbgWriter) = " .. activateDbgWriter)

print(" Parallelisation related parameters chosen:")
print("    distType   = " .. distributionType)

--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

-- create Instance of a Domain
print("Create Domain.")
dom = Domain()

-- load domain
print("Load Domain from File.")
LoadDomain(dom, gridName)

-- create Refiner
print("Create Refiner")
-- Create a refiner instance. This is a factory method
-- which automatically creates a parallel refiner if required.
refiner = GlobalDomainRefiner(dom)

-- Performing pre-refines
print("Perform (non parallel) pre-refinements of grid")
for i=1,numPreRefs do
	write( "PreRefinement step " .. i .. " ...")
	refiner:refine()
	print( " done.")
end

numProcs = GetNumProcesses()

-- Distribute the domain to all involved processes
-- Since only process 0 loaded the grid, it is the only one which has to
-- fill a partitionMap (but every process needs one and has to return his map
-- by calling 'DistributeDomain()', even if in this case the map is empty
-- for all processes but process 0).
if numProcs > 1 then
	util.DistributeDomain(dom, distributionType)	
end


--------------------------------------------------------------------------------
-- end of partitioning
--------------------------------------------------------------------------------

-- Perform post-refine
print("Refine Parallel Grid")
for i=numPreRefs+1,numRefs do
	write( "Refinement step " .. i .. " ...")
	refiner:refine()
	print( " done!")
end


-- Make sure, that the required subsets are present
requiredSubsets = {"Inner", "Boundary"}
if util.CheckSubsets(dom, requiredSubsets) == false then 
   print("Subsets missing. Aborting")
   exit()
end

print("NumProcs is " .. numProcs .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("ux", "Lagrange", 1)          
approxSpace:add_fct("uy", "Lagrange", 1)
approxSpace:add_fct("uz", "Lagrange", 1)    
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_local_dof_statistic(2)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-- lets order indices using Cuthill-McKee
OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Assembling
--------------------------------------------------------------------------------
print ("Setting up Assembling")

--------------------------------
-- User Data Functions
--------------------------------
function DirichletValue2d(x, y, t)
	return 0.0
end

function DirichletValue3d(x, y, z, t)
	return true, 0.0
end
	
--------------------------------------------------------------------------------
--  Setup FE Linear Element Discretization
--------------------------------------------------------------------------------

if dim == 3 then
	elemDisc = LinearElastoPlasticity("ux, uy, uz", "Inner") --Lame-equation
elseif dim == 2 then
	elemDisc = LinearElastoPlasticity("ux, uy", "Inner")
end

--elemDisc:set_elasticity_tensor_orthotropic(C11,C12,C13,C22,C23,C33,C44,C55,C66)
elemDisc:set_hooke_elasticity_tensor(11,8)
-- set volume force "f":
elemDisc:set_volume_forces(0.0) 
elemDisc:set_quad_order(2) -- this is a hack to use Prism, since no quad order > 2 available. Uncomment if not using Prism

dirichletBND = DirichletBoundary()

dirichletBND:add("DirichletValue"..dim.."d", "ux", "Boundary")
dirichletBND:add("DirichletValue"..dim.."d", "uy", "Boundary")
dirichletBND:add("DirichletValue"..dim.."d", "uz", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)






-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(50)
convCheck:set_minimum_defect(1e-15)
convCheck:set_reduction(1e-20)

print("done.")
-- create Linear Solver

-- create GMG solver --




-- Apply Solver

print ("Setting up Algebra Solver")

-- matrix and vectors
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
b2 = GridFunction(approxSpace)

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator


gs = GaussSeidel()
sgs = SymmetricGaussSeidel()

if true then
--[[	dbgWriter = GridFunctionDebugWriter(approxSpace)
	dbgWriter:set_base_dir(matrixWritePath)
	dbgWriter:set_vtk_output(false) ]]--
	transfer = StdTransfer(approxSpace)
	--transfer:set_debug(dbgWriter)
	
---------------------------------------------------------------
	-- Gemoetric Multi Grid
	
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_transfer(transfer)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0) -- 1
	gmg:set_base_solver(LU())
	gmg:set_smoother(sgs)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)

---------------------------------------------------------------

	amg = util.SetupRSAMGPreconditioner()
	
---------------------------------------------------------------
	
	gmgSolver = LinearSolver()
	
	gmgSolver:set_preconditioner(amg)
	--gmgSolver:set_preconditioner(Jacobi(0.5)) -- 3.62e0
	-- gmgSolver:set_preconditioner(Jacobi(0.6)) -- 5.86e+2 (rate 1.004 ! -> divergenz)
	--gmgSolver:set_preconditioner(GaussSeidel()) -- 8.790e-1
	-- gmgSolver:set_preconditioner(SymmetricGaussSeidel()) -- divergiert (rate 1.2e+2 !)
	gmgCGSolver = CG()
	gmgCGSolver:set_preconditioner(amg)
	
	-- create Convergence Check
	gmgConvCheck = ConvCheck()
	gmgConvCheck:set_maximum_steps(50)
	gmgConvCheck:set_minimum_defect(1e-24)
	gmgConvCheck:set_reduction(1e-24)
	
	gmgSolver:set_convergence_check(gmgConvCheck)
	gmgCGSolver:set_convergence_check(gmgConvCheck)
	
	baseSolver = LinearSolver()
	baseSolver:set_preconditioner(ILUT(0))
	baseSolver:set_convergence_check(gmgConvCheck)
	
	
---------------------------------------------------------------
	
	
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
	
	SaveMatrixForConnectionViewer(u, linOp, "Stiffness.mat")	
	
	
	b2 = GridFunction(approxSpace)
	b2:assign(b)
	
	iSeed = 4
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	domainDisc:adjust_solution(u)
	
	jac = Jacobi(0.5)
		
	if false then	
		--amg:set_debug_write_all_vectors(true)	
		gmgConvCheck:set_maximum_steps(40)
		print("Apply solver.")
		if ApplyLinearSolver(linOp, u, b, gmgSolver) == false then
			print("Could not apply linear solver.");
		end
	else		
			-- 2. apply solver
		print("Apply solver.")
		if ApplyLinearSolver(linOp, u, b, gmgCGSolver) == false then
			print("Could not apply linear solver.");
		end
		
		print("done");
		
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
		gmgConvCheck:set_maximum_steps(50)
		
		
		--gmg:set_debug(dbgWriter)
		
		
		if ApplyLinearSolver(linOp, u, b, gmgSolver) == false then
			print("Could not apply linear solver.");
		end
		SaveVectorDiffForConnectionViewer(u, uEx, "gmg_diff.vec")
		print("done");
		
		
		SaveVectorForConnectionViewer(u, "gmg_u.vec")
		--SaveVectorForConnectionViewer(u, A, "gmg_u.vec")	
		SaveVectorForConnectionViewer(b, "gmg_b.vec")
	end
end

WriteProfileData("famgLE3d.pdxml")