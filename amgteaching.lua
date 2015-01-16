----------------------------------------------------------
--
--   Lua - Script performing several test problems
--
--   Author: 	Arne Naegel
--   
-- 				based on laplace.lua by Andreas Vogel 
----------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("amgteaching_problems.lua")
--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:
dim = util.GetParamNumber("-dim", 2)

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	--gridName = "unit_square_tri_four_dirichlet_nodes.ugx"
	--gridName = "unit_square/unit_square_quads_8x8.ugx"
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
end

-- refinements:
numPreRefs = util.GetParamNumber("-numPreRefs", 2)
numRefs    = util.GetParamNumber("-numRefs",    5)

-- parallelisation related stuff
-- way the domain / the grid will be distributed to the processes:
distributionType = util.GetParam("-distType", "bisect") -- [grid2d | bisect | metis]

-- amount of output
verbosity = util.GetParamNumber("-verb", 0)	    -- set to 0 i.e. for time measurements,
						    -- >= 1 for writing matrix files etc.

activateDbgWriter = 0	  
activateDbgWriter = util.GetParamNumber("-dbgw", 0) -- set to 0 i.e. for time measurements,
						    -- >= 1 for debug output: this sets
						    -- 'fetiSolver:set_debug(dbgWriter)'

ilubeta = util.GetParamNumber("-ilu-beta", 0.0)

-- Display parameters (or defaults):
print(" General parameters chosen:")
print("    dim        = " .. dim)
print("    grid       = " .. gridName)
print("    numRefs    = " .. numRefs)
print("    numPreRefs = " .. numPreRefs)

print("    verb (verbosity)         = " .. verbosity)
print("    dbgw (activateDbgWriter) = " .. activateDbgWriter)

print(" Parallelisation related parameters chosen:")
print("    distType   = " .. distributionType)

print("    epsx = " .. epsx)
print("    epsy = " .. epsy)
print("    alpha = " .. alpha)
print("    ilu-beta = " .. ilubeta)
--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));

dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, {})
numProcs = NumProcs()

-- Distribute the domain to all involved processes
-- Since only process 0 loaded the grid, it is the only one which has to
-- fill a partitionMap (but every process needs one and has to return his map
-- by calling 'DistributeDomain()', even if in this case the map is empty
-- for all processes but process 0).
	
-- get subset handler
sh = dom:subset_handler()
if sh:num_subsets() ~= 2 then 
print("Domain must have 2 Subsets for this problem.")
exit()
end
sh:set_subset_name("Inner", 0)
sh:set_subset_name("DirichletBoundary", 1)
--sh:set_subset_name("NeumannBoundary", 2)

-- write grid to file for test purpose
if verbosity >= 1 then
	refinedGridOutName = "refined_grid_p" .. ProcRank() .. ".ugx"
	print("saving domain to " .. refinedGridOutName)
	SaveDomain(dom, refinedGridOutName)
	
	hierarchyOutName = "hierachy_p" .. ProcRank() .. ".ugx"
	print("saving hierachy to " .. hierarchyOutName)
	if SaveGridHierarchy(dom:grid(), hierarchyOutName) == false then
		print("Saving of domain to " .. hierarchyOutName .. " failed. Aborting.")
		    exit()
	end
end

print("NumProcs is " .. numProcs .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--  Order indices using Cuthill-McKee
-- OrderCuthillMcKee(approxSpace, true);


--Order indices lexicographically
OrderLex(approxSpace, "ul");

-- create algebraic Preconditioner
jac = Jacobi()
jac:set_damp(0.66)
gs = GaussSeidel()
bgs = BackwardGaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
--ilu:set_debug(dbgWriter)
ilu:set_beta(ilubeta)
ilut = ILUT()

------------------------------------------
-- solvers for alebraically smooth error
------------------------------------------
-- create convergence check (solver)
convCheckSmoother = ConvCheck()
convCheckSmoother:set_maximum_steps(5)
convCheckSmoother:set_reduction(1e-12)
convCheckSmoother:set_minimum_defect(1e-14)
convCheckSmoother:set_verbose(false)
convCheckSmoother:set_supress_unsuccessful(true)


-- create Jacobi Solver
jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheckSmoother)

-- create Gauss-Seidel Solver
gsSolver = LinearSolver()
gsSolver:set_preconditioner(gs)
gsSolver:set_convergence_check(convCheckSmoother)

-- create SGS Solver
sgsSolver = LinearSolver()
sgsSolver:set_preconditioner(sgs)
sgsSolver:set_convergence_check(convCheckSmoother)

-- create ILU Solver
iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilu)
iluSolver:set_convergence_check(convCheckSmoother)


-- select smoother/solver

smoother={}
--smoother.solver = jacSolver
--smoother.name = "Jacobi"
--smoother.pre = jac
--smoother.post = jac
smoother.solver = gsSolver
smoother.name = "GS"
smoother.pre = gs
smoother.post = bgs
--smoother.name = "SGS"
--smoother.solver = sgsSolver
--smoother.pre = sgs
--smoother.post = sgs
--smoother.solver = iluSolver
--smoother.name = "ILU"
--smoother.pre = ilu
--smoother.post = ilu

-- create (exact) base solver
-----------------

baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(1000)
baseConvCheck:set_minimum_defect(1e-11)
baseConvCheck:set_reduction(1e-30)
baseConvCheck:set_verbose(false)
base = LinearSolver()
base:set_convergence_check(baseConvCheck)
base:set_preconditioner(ilu)

base = LU()



-- create GMG ---
-----------------

-- Geometric Multi Grid
gmg = GeometricMultiGrid(approxSpace)

gmg:set_base_level(0)
gmg:set_base_solver(base)
gmg:set_gathered_base_solver_if_ambiguous(false)
gmg:set_smoother(smoother.pre)
gmg:set_cycle_type(1)
gmg:set_num_presmooth(1)
gmg:set_num_postsmooth(1)
if activateDbgWriter >= 1 then
gmg:set_debug(dbgWriter)
end

-- create AMG ---
-----------------
if false then
amg = RSAMGPreconditioner()
amg:set_nu1(2)
amg:set_nu2(2)
amg:set_gamma(1)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)
amg:set_base_solver(base)
--amg:set_debug(u)
end




------------------------------------------
-- real solvers
------------------------------------------
-- create convergence check (solver)
convCheckSolver = ConvCheck()
convCheckSolver:set_maximum_steps(50)
convCheckSolver:set_minimum_defect(1e-14)
convCheckSolver:set_reduction(1e-12)
convCheckSolver:set_verbose(false)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(convCheckSolver)

-- create CG Solver
cgSolver = CG()
cgSolver:set_preconditioner(ilu)
cgSolver:set_convergence_check(convCheckSolver)

-- create BiCGStab Solver
bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(jac)
bicgstabSolver:set_convergence_check(convCheckSolver)
solver = linSolver
--------------------------------------------------------------------------------
--  Setup User Functions
--------------------------------------------------------------------------------



--------------------------------------------------------------------------------
--  Run problems (with several configurations)
--------------------------------------------------------------------------------	

-- for each of the problems
for myname,myproblem in pairs(amgTeaching.problem) do
	
	
	
	elemDisc = amgTeaching.CreateProblemDependentElemDisc(myproblem, dim)
	print ("\n---------------------------------------------------------------------------------")
	print("Running "..myproblem.desc.." ["..myname.."]")
	print ("---------------------------------------------------------------------------------")
	--util.PrintTable(myproblem)
	
	
	-- for each configuration
	--print(myproblem.configs)
	for key,config in pairs(myproblem.configs) do
		print ("----- Running parameter set "..key..": -----")
		print ("  epsx="..config.epsx..", epsy="..config.epsy)
		epsx = config.epsx
		epsy = config.epsy
		
		domainDisc = amgTeaching.SetUpProblem(approxSpace, myproblem, elemDisc)
		A, x, b = amgTeaching.AssembleAxb(approxSpace, domainDisc)

		-- update GMG ---
		-----------------
		gmg:set_discretization(domainDisc)	
		
		--------------------------------------------------------------------------------
		--  Apply solver - using method defined in 'operator_util.h',
		--  to get separate profiling for assemble and solve
		--------------------------------------------------------------------------------
		
		-- 1.b write matrix for test purpose
		if verbosity >= 1 then
			SaveMatrixForConnectionViewer(u, A, myname.."P"..key.."Stiffness.mat")
			SaveVectorForConnectionViewer(b, myname.."P"..key.."Rhs.vec")
		end
		
		-- 2. apply solver
		--print("Apply solver.")
		ApplyLinearSolver(A, u, b, solver)
		
		print(" step: "..convCheckSolver:step())
		print(" avg_rate: "..convCheckSolver:avg_rate())
		print(" last reduction: "..convCheckSmoother:previous_defect() / convCheckSmoother:defect())
		
		if verbosity >= 1 then
			WriteGridFunctionToVTK(u, myname.."P"..key.."Solve"..smoother.name)
		end
	
		
		--------------------------------------------------------------------------------
		--  Apply smoother - this generates algebraically smooth errors
		--------------------------------------------------------------------------------

		-- set u to random (except dirichlet point), and b = 0		
		u:set_random(-1.0, 1.0)
		domainDisc:adjust_solution(u)
		b:set(0.0)
		
		-- 2. apply solver
		ApplyLinearSolver(A, u, b, smoother.solver)
		
		
		if verbosity >= 1 then
			WriteGridFunctionToVTK(u, myname.."P"..key.."Smooth"..smoother.name)
		end
	
	end -- for each configuration
end -- for each problem



--------------------------------------------------------------------------------
--  Print Profiling
--------------------------------------------------------------------------------
