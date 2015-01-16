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
ug_load_script("util/mac_util.lua")
--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------

-- command line parameters
--------------------------------------------------------------------------------------------------
--matrixWritePath=os.getenv("HOME").."/matrices/"

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

iTestvectorSmooths = util.GetParamNumber("-iTestvectorSmooths", 0, "nr of testvector smooth for construction of P")
bUseJacobi = util.HasParamOption("-bUseJacobi", "if true, use jacobi as smoother. otherwise use gs")
iSmooths = util.GetParamNumber("-iSmooths", 2, "nr of pre/postsmooths")

-- Several definitions which can be changed by command line parameters

dim = 3 -- dim is always 3

-- CEMBALO PARAMETERS --

--thePath=os.getenv("HOME").."/matrices/"	
thePath = util.createMacRAMDisk("ug4RAMDisk", 512).."/RAMmatrices/"
CreateDirectory(thePath) 
matrixWritePath=thePath

-- grid
gridName = util.GetParam("-grid", "cembalo/cembalo_3d_test2.ugx")
gridName = util.GetParam("-grid", "cembalo/rechteck_100_hoch.ugx")
gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
gridName = util.GetParam("-grid", "cembalo/cube_one_bnd.ugx")

gridName = util.GetParam("-grid", "cembalo/cembalo_3d_test2.ugx")

gridName = util.GetParam("-grid", "cembalo/cube_one_bnd.ugx")
gridName = util.GetParam("-grid", "amgtests/cube_one_bnd.ugx")

-- rho = density
rho = util.GetParamNumber("-rho", 0.00050, "density")

--SetDebugLevel(debugID.AMG.Coarsening, 4)
--SetDebugLevel(debugID.AMG.FAMG.Interpolation, 4)

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
InitUG(dim, AlgebraType("CPU", 3));

-- create Instance of a Domain
print("Create Domain.")

neededSubsets = {"Inner", "Boundary"}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

stretchX = util.GetParamNumber("-stretchX", 1.0)
stretchY = util.GetParamNumber("-stretchY", 1.0)
stretchZ = util.GetParamNumber("-stretchZ", 1.0)
print("scaling domain by ("..stretchX..", "..stretchY..", "..stretchZ..")") 
ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, stretchZ))



print("NumProcs is " .. NumProcs() .. ", numPreRefs = " .. numPreRefs .. ", numRefs = " .. numRefs .. ", grid = " .. gridName)

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

matLaw = HookeLaw()
matLaw:set_hooke_elasticity_tensor(11,8)
--matLaw:set_elasticity_tensor_orthotropic(C11,C12,C13,C22,C23,C33,C44,C55,C66)

if dim == 3 then
	elemDisc = SmallStrainMechanics("ux, uy, uz", "Inner") --Lame-equation
elseif dim == 2 then
	elemDisc = SmallStrainMechanics("ux, uy", "Inner")
end
elemDisc:set_material_law(matLaw)
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




	amg = FAMGPreconditioner()
	amg = util.SetupFAMGPreconditioner(SuperLU(), presmoother, postsmoother)
	amg:set_write_f_values(true)
	amg:set_verbose(1)
	
	amg:set_XC_LE_griebel(true)
	
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
	

	amg:set_max_nodes_for_base(maxBase)
	amg:set_max_levels(maxLevels)
	
	-- amg:set_debug_level_calculate_parent_pairs(2)
	-- amg:set_debug_level_precalculate_coarsening(4)
	-- amg:set_debug_level_external_coarsening(4)
	
	jac2 = Jacobi()
	jac2:set_damp(0.4)
	amg:set_testvector_smooths(0)
	amg:set_testvector_smoother(jac2)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
	
	if false then
		amg:set_testvector_from_matrix_rows()
	else
		-- amg:set_testvectors_LE(1e4)
		-- amg:set_testvectors_LE(1e7) -- gg
		--amg:set_testvectors_LE(1e3)
		bWriteMat = true		
	end	
	
	
	amg:set_delta(famgDelta)
	 
	amg:set_theta(0.8)
	if util.HasParamOption("-tvLin") then 
		amg:set_testvectors(TestvectorsLinear3d())
	elseif util.HasParamOption("-tvLE") then
		amg:set_testvectors(TestvectorsLinearElasticity3d())
	else
		amg:set_testvectors(TestvectorsConstant())
	end
	
	amg:set_one_init(true)
	
	vectorWriter = GridFunctionPositionProvider()
	vectorWriter:set_reference_grid_function(u)
	amg:set_position_provider(vectorWriter)

	amg:set_write_testvectors(true)
	amg:set_write_f_values(true)
	amg:set_matrix_write_path(matrixWritePath)


A = MatrixOperator()
-- init Operator and set dirichlet values in start iterate
print("Assemble A and b...")
tBefore = os.clock()
domainDisc:assemble_linear(A, b2)
tAssemble = os.clock()-tBefore
print ("done")

print(amg:config_string())

print("Init solver for operator.")
amgInit = LinearSolver()
amgInit:set_preconditioner(amg)
amgInit:init(A)	

print(amg:level_information_str())

results = {}


iSeed = 2

function solve_problem(linSolver, name)
	tBefore = os.clock()
	linSolver:init(A)
	srand(iSeed)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	
	linSolver:apply_return_defect(u,b)	
	tSolve = os.clock()-tBefore
	res  ={}
	res.lastReduction = convCheck:defect()/convCheck:previous_defect()
	res.avgReduction = convCheck:avg_rate()
	res.steps = convCheck:step()
	res.time = tSolve
	table.insert(results, {name, res})
	print(name.." took "..tSolve.." s, \tavgReduction: "..convCheck:avg_rate()..", \tlastReduction: "
	..res.avgReduction..", \tsteps = "..convCheck:step())
end
function solve_lin(linSolver, amg, name)
	linSolver:set_convergence_check(convCheck)	
	linSolver:set_preconditioner(amg)
	solve_problem(linSolver, name)
end

function solve_amg(smoother, numSmooths, name)

	amg:set_presmoother(smoother)
	amg:set_postsmoother(smoother)
	amg:set_num_presmooth(numSmooths)
	amg:set_num_postsmooth(numSmooths)
	
	linSolver = LinearSolver()
	linSolver:set_preconditioner(amg)
	cgSolver = CG()
	cgSolver:set_preconditioner(amg)
	
	
	amg:init_smoothers()
	
	print(name..": "..numSmooths.."x "..smoother:config_string())
	
	
	--ScaleDomain(dom, MakeVec(0,0,0), MakeVec(1/stretchX, 1/stretchY, 1.0))
	--SaveVectorForConnectionViewer(u, "amg_u_linear.vec")
	--ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))
		
	print("CG, "..name)
	solve_lin(CG(), amg, "CG("..name..")")
	print("LinearSolver, "..name)
	solve_lin(LinearSolver(), amg, "LS("..name..")")	
	

end

convCheck:set_verbose(false)
solve_amg(SymmetricGaussSeidel(), 2, "GSx2")

amg:check_vector(TestvectorsLinear3d(), "LIN")
amg:check_vector(TestvectorsLinearElasticity3d(), "LE")



solve_amg(SymmetricGaussSeidel(), 4, "GSx4")
solve_amg(SymmetricGaussSeidel(), 10, "GSx10")
solve_amg(SymmetricGaussSeidel(), 20, "GSx20")


if util.HasParamOption("-check") then
amg:set_fsmoothing(false)
	amg:set_num_presmooth(2)
	amg:set_num_postsmooth(2)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	linSolver = LinearSolver()
	linSolver:set_preconditioner(SymmetricBlockGaussSeidel())
	linSolver:init(A)		
	convCheck:set_maximum_steps(0)
	linSolver:set_convergence_check(convCheck)	
	linSolver:apply_return_defect(u,b)
	--amg:check(u, b)
	amg:check()
end



WriteProfileData("famgLE.pdxml")
SaveVectorForConnectionViewer(u, "Solution.vec")
