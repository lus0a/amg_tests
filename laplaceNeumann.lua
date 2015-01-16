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
--gridName = util.GetParam("-grid", "unit_square_quads_8x8_four_dirichlet_nodes.ugx")
gridName = util.GetParam("-grid", "unit_square/neumannGrid.ugx")
-- gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_dirichlet_bottom.ugx")


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
approxSpace:add_fct("c", "Lagrange", 1)
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


function ourNeumannBnd2d(x, y, t)
   return true, -1
end

-- dirichlet setup
dirichlet = LuaCondUserNumber("ourDirichletBnd"..dim.."d")
neumann = LuaCondUserNumber("ourNeumannBnd"..dim.."d")

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

	diffusionMatrix = ConstUserMatrix(1.0)
	dirichlet = LuaCondUserNumber("ourDirichletBnd2d")
		
	function Diffusion2d(x, y, t)
		return	1, 0, 
				0, 1
	end
	elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
	elemDisc:set_diffusion("Diffusion2d")
	
	dirichletBND = DirichletBoundary()
	dirichletBND:add(dirichlet, "c", "Boundary")
	
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(elemDisc)
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
	amg = util.SetupFAMGPreconditioner(base, presmoother, postsmoother)
	
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
	amg:set_delta(0.127)
	amg:set_delta(0.11) -- g
	-- amg:set_delta(1e7)
	-- amg:set_delta(1e5)		
	amg:set_theta(0.9)
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
	
	if false then
		amg:set_testvector_from_matrix_rows()
	else
		-- amg:set_testvectors_LE(1e4)
		-- amg:set_testvectors_LE(1e7) -- gg
		amg:set_testvectors_LE(1)
		bWriteMat = true		
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

		


-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(50)
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
	transfer = StdTransfer(approxSpace)
	transfer:set_debug(dbgWriter)
	
	
	-- Gemoetric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_transfer(transfer)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0) -- 1
	gmg:set_base_solver(base)
	gmg:set_smoother(gs)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(5)
	gmg:set_num_postsmooth(5)
	gmgSolver = LinearSolver()
	gmgSolver:set_preconditioner(gmg)
	gmgCGSolver = CG()
	gmgCGSolver:set_preconditioner(gmg)
	
	-- create Convergence Check
	gmgConvCheck = ConvCheck()
	gmgConvCheck:set_maximum_steps(200)
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
	gmgConvCheck:set_maximum_steps(20)
	
	
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

----
print("\n")		
util.PrintArguments()
util.PrintIgnoredArguments()