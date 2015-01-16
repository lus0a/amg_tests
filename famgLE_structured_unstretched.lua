--------------------------------------------------------------------------------
--[[!
-- \file apps/amg/famgLE_structured_unstretched.lua
--
-- test ob FAMG LE auf einheitsquadrat, mit eingespannt an einer Seite
--
]]--
--------------------------------------------------------------------------------

-- ProfileLUA(true)
--SetOutputProfileStats(true)
InitSignals()
ug_load_script("ug_util.lua")
ug_load_script("util.lua")
ug_load_script("solver_util/setup_famg.lua")
ug_load_script("util/mac_util.lua")
ug_load_script("setup.lua")
-- command line parameters

--InitSignals()	
--------------------------------------------------------------------------------------------------

paramOutputProcess = util.GetParamNumber("-outproc", 0, "Output processor")
numRefs = util.GetParamNumber("-numRefs", 1, "number of refinements")
numPreRefs = util.GetParamNumber("-numPreRefs", -1, "number of refinements before ")
maxBase = util.GetParamNumber("-maxBase", 1000, "maximum number of nodes before base solver is used")
maxLevels = util.GetParamNumber("-maxLevels", 2, "maximum number of AMG levels")

bWriteStats = util.HasParamOption("-bWriteStats", "if true, write stats output")
bWriteMat =util.HasParamOption("-writeMatrices", "if true, write matrix output")

bAggressiveCoarsening = util.HasParamOption("-AC", "if true, use aggressive coarsening")
famgDelta = util.GetParamNumber("-delta", 0.2, "Interpolation quality F may not be worse than this (F < m_delta)")
famgTheta = util.GetParamNumber("-theta", 0.9, "with multiple parents, discard pairs with m_theta * F > min F.")

logtofile = util.GetParam("-logtofile", nil, "the file to write logs to")
bWCycle = util.HasParamOption("-Wcycle", "if true, use W-Cycle, otherwise: V-Cycle.")

iTestvectorSmooths = util.GetParamNumber("-iTestvectorSmooths", 0, "nr of testvector smooth for construction of P")
bUseJacobi = util.HasParamOption("-bUseJacobi", "if true, use jacobi as smoother. otherwise use gs")
iSmooths = util.GetParamNumber("-iSmooths", 2, "nr of pre/postsmooths")


if util.HasParamOption("-special") then
	gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_4bnd2.ugx")  addNumRefs = 0
elseif util.HasParamOption("-33k") then
	BoundaryName = "Bottom"
	gridName = util.GetParam("-grid", "amgtests/unstructured_33k.ugx") --- numrefs+1 - d 0.03
elseif util.HasParamOption("-3k") then
	BoundaryName = "Bottom"
	gridName = util.GetParam("-grid", "amgtests/unstructured_3k.ugx")
elseif util.HasParamOption("-cembalo") then
	gridName = util.GetParam("-grid", "cembalo/cembalo_platte_2d_sci_1.ugx")
	numRefs =0
	numPreRefs =0
	BoundaryName = "Boundary"
else
	BoundaryName = "Bottom"
	gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_4bnd.ugx")
end

WriteProfileData("famgLE.pdxml")
if false then
	gridName = util.GetParam("-grid", "amgtests/beam2d.ugx")
	gridName = util.GetParam("-grid", "beamlong2.ugx")
	BoundaryName = "Boundary"
end

--------------------------------------------------------------------------------------------------

function check()
	print("-----------")
	amg:set_fsmoothing(false)
	amg:set_num_presmooth(5)
	amg:set_num_postsmooth(5)
	print(amg:config_string())
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
	print("-----------")
end


stringResults=""
function solve_problem(linSolver, name)
	tBefore = os.clock()
	linSolver:init(A)
	srand(iSeed)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	SaveVectorForConnectionViewer(u, "u.vec")
	SaveVectorForConnectionViewer(b, "b.vec")
	linSolver:apply_return_defect(u,b)	
	tSolve = os.clock()-tBefore
	res  ={}
	res.lastReduction = convCheck:defect()/convCheck:previous_defect()
	res.avgReduction = convCheck:avg_rate()
	res.steps = convCheck:step()
	res.time = tSolve
	table.insert(results, {name, res})
	print(name.." took "..tSolve.." s, \tavgReduction: "..convCheck:avg_rate()..", \tlastReduction: "
	..res.lastReduction..", \tsteps = "..convCheck:step().."\ttSolve = "..tSolve)
	stringHeader="name\ttook\tavgReduction\tlastReduction\tsteps\tsolve"
	stringResults=stringResults..name.."\t"..tSolve.."\t"..convCheck:avg_rate().."\t"
	..res.lastReduction.."\t"..convCheck:step().."\t"..tSolve.."\n"
end
function solve_lin(linSolver, amg, name)
	linSolver:set_convergence_check(convCheck)	
	linSolver:set_preconditioner(amg)
	
	solve_problem(linSolver, name)
	
end
function printStringResults()
	print(stringHeader)
	print(stringResults)
end

function solve_amg(smoother, numSmooths, name)


	amg:set_presmoother(smoother)
	amg:set_postsmoother(smoother)
	amg:set_num_presmooth(numSmooths)
	amg:set_num_postsmooth(numSmooths)
	
	linSolver = LinearSolver()
	linSolver:set_preconditioner(amg)
	cgSolver = BiCGStab()
	cgSolver:set_preconditioner(amg)
	
	
	amg:init_smoothers()
	
	print(name..": "..numSmooths.."x "..smoother:config_string())
	
	
	--ScaleDomain(dom, MakeVec(0,0,0), MakeVec(1/stretchX, 1/stretchY, 1.0))
	--SaveVectorForConnectionViewer(u, "amg_u_linear.vec")
	--ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))
		
	print("CG, "..name)
	--solve_lin(BiCGStab(), amg, "BiCGStab("..name..")")
	solve_lin(LinearSolver(), amg, "LS("..name..")")
	--print("LinearSolver, "..name)
	--solve_lin(LinearSolver(), amg, "LS("..name..")")	
	

end


--------------------------------------------------------------------------------------------------
util.CheckAndPrintHelp(
"Script for FAMG on Linear Elasticity\nMartin Rupp"
)

util.bFailOnMultipleParameters= true

--------------------------------------------------------------------------------------------------
-- set up output
if bWriteStats then
	outdir=os.getenv("HOME").."/results/famgLE/"
	commandLineFilename = util.GetUniqueFilenameFromCommandLine()
	thePath = outdir.."famgLE"..commandLineFilename.."/"
	CreateDirectory(thePath)
	
	GetLogAssistant():enable_file_output(false, "")
	GetLogAssistant():enable_file_output(true, thePath.."log.txt")
else
	--thePath=os.getenv("HOME").."/matrices/"	
	thePath = util.createMacRAMDisk("ug4RAMDisk", 768).."/RAMmatrices/"
	CreateDirectory(thePath) 
end

matrixWritePath=thePath

a = GetLogAssistant()
a:set_output_process(-1)
if paramOutputProcess ~= -1 and ProcRank() ~= paramOutputProcess then a:enable_terminal_output(false) end
a:enable_file_output(true, thePath.."log_p"..ProcRank()..".txt")
print("Writing matrices to "..matrixWritePath)
if bWriteMat and outdir == "" then
	outdir = matrixWritePath
end

--------------------------------------------------------------------------------------------------

dim = 2

-- choose algebra
bLinearElasticity = true
if util.HasParamOption("-laplace") then
	bLinearElasticity = false
end

if bLinearElasticity then
InitUG(dim, AlgebraType("CPU", 2))
else
InitUG(dim, AlgebraType("CPU", 1))
end

-- Create, Load, Refine and Distribute Domain
neededSubsets = {"Inner", BoundaryName}
dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

stretchX = util.GetParamNumber("-stretchX", 1.0)
stretchY = util.GetParamNumber("-stretchY", 1.0)
print("scaling domain by ("..stretchX..", "..stretchY..")") 
ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))

d = util.GetParamNumber("-d", 0.1)
if d > 0.0 then
	d = math.pow(d, numRefs+1)
	print("randomize ="..d)
	RandomizeDomain(dom, d, d, d)
end

print("stretchX = "..stretchX..", stretchY = "..stretchY)

-- create Approximation Space
print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom)
-------------------------------------------
--  Setup User Functions
-------------------------------------------
print ("Setting up Assembling")
-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------

if bLinearElasticity then

	approxSpace:add_fct("ux", "Lagrange", 1)          
	approxSpace:add_fct("uy", "Lagrange", 1)
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()
	
	elemDisc = SmallStrainMechanics("ux, uy", "Inner")
	
	matLaw = HookeLaw()
	-- set Elasticity-tensor "C" with material constants: 
	-- first parameter: lambda ("1.lame constant")
	-- second parameter: mu ("2.lame constant")
	if util.HasParamOption("-bSmall") then
		--matLaw:set_hooke_elasticity_tensor(0.11074382,0.08019380)
	else
		--matLaw:set_hooke_elasticity_tensor(11,8)
	end
	
	nu = util.GetParamNumber("-nu", 0.4)
	matLaw:set_hooke_elasticity_tensor_E_nu(20.6316, nu)
	elemDisc:set_material_law(matLaw)
else
	approxSpace:add_fct("c", "Lagrange", 1)          
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()
	
	
	elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
	elemDisc:set_diffusion(1)
	elemDisc:set_source(0)

end
-----------------------------------------------------------------
--  Setup Dirichlet Boundary
-----------------------------------------------------------------

function DirichletValue2d(x, y, t)
	return 1
end

dirichletBND = DirichletBoundary()
function setBoundary(boundaryName)
	if bLinearElasticity then
		dirichletBND:add("DirichletValue"..dim.."d", "ux", boundaryName)
		dirichletBND:add("DirichletValue"..dim.."d", "uy", boundaryName)
	else
		dirichletBND:add("DirichletValue"..dim.."d", "c", boundaryName)
	end
end

bnd = util.GetParamNumber("-bnd", 1)
setBoundary(BoundaryName)
if bnd > 1 then
setBoundary("Top")
if bnd > 2 then
setBoundary("Left")
if bnd > 3 then
setBoundary("Right")
end end end


-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)
-- domainDisc:add(neumannDisc)

-- ? not necessary
OrderCuthillMcKee(approxSpace, false);

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

SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
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

--postsmoother = BlockGaussSeidel(1)
postsmoother = GaussSeidel()

presmoother = postsmoother




base = LU()
base:set_minimum_for_sparse(0)

------- Create FAMG Solver ---------
amg = FAMGPreconditioner()
--amg = util.SetupFAMGPreconditioner(base, presmoother, postsmoother)
amg:set_num_presmooth(iSmooths)
amg:set_num_postsmooth(iSmooths)
amg:set_base_solver(SuperLU())
amg:set_min_nodes_on_one_processor(1)
amg:use_checks(true)

--- ATTENTION
amg:set_xc_calc_interpolation_routine(2)
--- ATTENTION

--SetDebugLevel(debugID.AMG.Agglomeration, 5)

--SetDebugLevel(debugID.AMG.Coarsening, 4)
--SetDebugLevel(debugID.AMG.FAMG.Interpolation, 6)


--amg:add_watchnode(0,1006)
--amg:add_watchnode(0, 883)
amg:set_indirect_interpolation(false)
if bWCycle then
	amg:set_cycle_type(2)
else
	amg:set_cycle_type(1)
end

amg:set_write_f_values(true)



	
amg:set_max_nodes_for_base(maxBase)
amg:set_max_levels(maxLevels)

amg:set_strong_connection_external(0.3)
--amg:set_strong_connection_submask(1)

amg:set_testvector_smooths(iTestvectorSmooths)
amg:set_testvector_smoother(Jacobi(0.3))
--amg:set_testvector_smoother(BlockGaussSeidel(2))
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)

if util.HasParamOption("-XC") then
amg:set_external_coarsening(true)
end



if util.HasParamOption("-tvLin") then
tvDesc = "tvLin"
amg:set_testvectors(TestvectorsLinear2d())
elseif util.HasParamOption("-tvL") then
tvDesc = "tvLComp"
amg:set_testvectors(TestvectorsLinearComponent())
elseif util.HasParamOption("-tvLE") then
tvDesc = "tvLE"
amg:set_testvectors(TestvectorsLinearElasticity2d())
elseif util.HasParamOption("-tvAuto") then
tvDesc = "tvAuto"
--amg:set_testvectors(TestvectorsAuto(10, Jacobi(0.3), 50))
amg:set_testvectors(TestvectorsAuto(30, 2, SymmetricGaussSeidel(), 300))
elseif util.HasParamOption("-tvCWobbled") then
amg:set_testvectors(TestvectorsConstantWobbled(0.1))
elseif util.HasParamOption("-tvOR") then
tvDesc = "tvOR"
amg:set_testvectors(TestvectorsOnlyRotation2d())
else
tvDesc = "tvConst"
amg:set_testvectors(TestvectorsConstant())
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
if util.HasParamOption("-XCgriebel") then
	amg:set_XC_LE_griebel(true)
end

if util.HasParamOption("-XCnew") then
	amg:set_XC_new(true)
end


amg:set_one_init(true)
amg:set_injection_restriction( util.HasParamOption("-bUseInjection") )

if bWriteMat then
	amg:set_write_testvectors(true)
	amg:set_write_f_values(true)
	amg:set_matrix_write_path(matrixWritePath)
end	

amg:set_deficit_solve_tolerance(1e-7)
amg:set_deficit_solve_tolerance(1e-3)
--amg:set_deficit_solve_tolerance(1e-7)
--amg:set_deficit_solve_tolerance(1e-2)



amg:set_position_provider( GridFunctionPositionProvider(u))

-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(40)
convCheck:set_minimum_defect(1e-12)
convCheck:set_reduction(1e-20)
convCheck:set_verbose(true)	


if util.HasParamOption("-bYcycle") then
	convCheckY = ConvCheck()
	convCheckY:set_maximum_steps(4)
	convCheckY:set_minimum_defect(1e-30)
	convCheckY:set_reduction(0.1)
	convCheckY:set_verbose(false)
	
	--amg:set_Y_check(convCheckY)
	amg:set_Y_cycle(40, 1e-30, 0.01)	
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

-- Geometric Multi Grid
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)
gmg:set_base_solver(LU())
gmg:set_smoother(BlockGaussSeidel(1));
gmg:set_num_presmooth(5)
gmg:set_num_postsmooth(5)
gmg:set_cycle_type(1)

if false then
	linSolver = LinearSolver()
	linSolver:set_preconditioner(SymmetricBlockGaussSeidel())
	tBefore = os.clock()
	linSolver:init(A)
	srand(iSeed)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	linSolver:apply_return_defect(u,b)	
	SaveVectorForConnectionViewer(u, "u.vec")
	SaveVectorForConnectionViewer(b, "b.vec")
end
	
	



results = {}

if util.HasParamOption("-adaptive") then
	amg:set_num_presmooth(iSmooths)
	amg:set_num_postsmooth(iSmooths)
	print(amg:config_string())
	amg:set_one_init(false)
	
	--tva = TestvectorsAuto(10, 1, SymmetricGaussSeidel(), 10)
	tva = TestvectorsConstant()
	amg:adaptive_filtering(A, u, ConvCheck(14,1e-12, 1e-8, true), 0.4, tva)
	--amg:adaptive_filtering(A, u, convCheck, 0.4, TestvectorsConstant())
	do return end
end


-- choose precond
print("Init solver for operator.")
amgInit = LinearSolver()
amgInit:set_preconditioner(amg)
amgInit:init(A)	

print(amg:level_information_str())

print(amg:config_string())
--solve_amg(SymmetricGaussSeidel(), 2, "SGSx2")


if( NumProcs() == 1 and bLinearElasticity) then
amg:check_vector(TestvectorsLinear2d(), "LIN")
amg:check_vector(TestvectorsLinearElasticity2d(), "LE")
end

solve_amg(SymmetricGaussSeidel(), 1, "SGSx01")
solve_amg(SymmetricGaussSeidel(), 2, "SGSx02")
solve_amg(SymmetricGaussSeidel(), 4, "SGSx04")
solve_amg(SymmetricGaussSeidel(), 8, "SGSx08")
solve_amg(SymmetricGaussSeidel(), 16, "SGSx16")
solve_amg(SymmetricGaussSeidel(), 32, "SGSx32")
solve_amg(SymmetricGaussSeidel(), 64, "SGSx64")
solve_amg(SymmetricGaussSeidel(), 128, "SGSx128")

if false then
solve_amg(Jacobi(0.6), 1, "Jacx01")
solve_amg(Jacobi(0.6), 2, "Jacx02")
solve_amg(Jacobi(0.6), 4, "Jacx04")
solve_amg(Jacobi(0.6), 8, "Jacx08")
solve_amg(Jacobi(0.6), 16, "Jacx16")
solve_amg(Jacobi(0.6), 32, "Jacx32")
end

--solve_amg(SymmetricGaussSeidel(), 20, "SGSx20")
printStringResults()
if false then
	do return end
	--amg:set_fsmoothing(false)
	
	solve_amg(SymmetricGaussSeidel(), 50, "SGSx50")
	
	
	solve_amg(SymmetricBlockGaussSeidel(1), 10, "SBGS(1)x10")
	--check()
	do return end
	solve_amg(ILU(), 2, "ILU")
	
	
	
	if util.HasParamOption("-check") then
		check()
	end
end


WriteProfileData("famgLE.pdxml")
SaveVectorForConnectionViewer(u, "Solution.vec")

--if true then return end
--solve_amg(SymmetricGaussSeidel(), 50, "GSx50")


--solve_problem(BiCGStab(ILU()), "ILLU")
--solve_problem(SuperLU(), "SuperLU")

if false then
	hypre = HYPRE()
	dbgWriter = GridFunctionDebugWriter(approxSpace)
	dbgWriter:set_vtk_output(true)
	hypre:set_debug(dbgWriter)
	solve_problem(hypre, "HYPRE")
end	

SaveVectorForConnectionViewer(u, "Solution.vec")

if ProcRank() == 0 then
	
	amgLevels = amg:get_used_levels()
	amgInfo={}
	for i=0,amgLevels-1 do
		amgInfo[i] = amg:get_level_information(i):get_nr_of_nodes()
	end
	
	stats = {	
	{ "gridName", gridName},
	{ "numProcs", NumProcs() },
	{ "numRefs", numRefs},
	{ "bnd", bnd},
	{ "stretchX", stretchX},
	{ "nu", nu},
	{ "TV", tvDesc},	
	{ "delta", famgDelta},
	{ "opComp", amg:get_operator_complexity()},
	{ "gridComp", amg:get_grid_complexity()},
	
	}
	for i,r in ipairs(results) do
		name = r[1]
		res = r[2]
		print(name)		
		table.insert(stats, {"name", name})
		table.insert(stats, {"lastReduction", res.lastReduction})
		table.insert(stats, {"avg", res.avgReduction})
		table.insert(stats, {"steps", res.steps})
		table.insert(stats, {"time", res.time})
	end
	flattenTable(getAMGLevelTable(amg), "level", stats)
	--util.PrintTable(stats)
	--util.printStats(stats)
	--if bWriteStats  then	
		util.writeFileStats(stats, "/Users/mrupp/results/mystats.txt")
	--end
end


--do return end

-- fix SGS!
--solve_amg(SymmetricGaussSeidel(), 2, "lin_SGS_2")
--solve_amg(SymmetricGaussSeidel(), 4, "lin_SGS_4")
--solve_amg(SymmetricGaussSeidel(), 10, "lin_SGS_10")
do return end
solve_amg(BlockGaussSeidel(1), 1, "lin_BGS1_1")
solve_amg(BlockGaussSeidel(1), 2, "lin_BGS1_2")
solve_amg(BlockGaussSeidel(2), 2, "lin_BGS2_2")

do return end
solve_amg(BlockGaussSeidel(3), 2, "lin_BGS3_2")
solve_amg(SymmetricBlockGaussSeidel(1), 1, "lin_SBGS1_1")
solve_amg(SymmetricBlockGaussSeidel(1), 2, "lin_SBGS1_2")
solve_amg(BlockGaussSeidelIterative(1, 10), 1, "lin_BGSIt1_1")
solve_amg(BlockGaussSeidelIterative(1, 10), 2, "lin_BGSIt1_2")
solve_amg(BlockGaussSeidelIterative(2, 10), 2, "lin_BGSIt2_2")
solve_amg(BlockGaussSeidelIterative(3, 10), 2, "lin_BGSIt3_2")
solve_amg(SymmetricBlockGaussSeidelIterative(1, 10), 1, "lin_SBGSIt1_1")
solve_amg(SparseBlockGaussSeidel(3), 1, "lin_SBGSIt1_1")
solve_amg(SparseBlockGaussSeidel2(3), 1, "lin_SBGSIt1_1")


do return end

if false then
solve_amg(ILUTScalar(10000), 2, "lin_ILUT_100")
--solve_amg(ILUTScalar(10), 2, "lin_ILUT_10")
---solve_amg(ILUTScalar(1), 2, "lin_ILUT_1")
solve_amg(ILUTScalar(0.1), 2, "lin_ILUT_0.1")
solve_amg(ILUTScalar(1e-2), 2, "lin_ILUT_1e-2") ----
solve_amg(ILUTScalar(1e-3), 2, "lin_ILUT_1e-3") ----
solve_amg(ILUTScalar(1e-4), 2, "lin_ILUT_1e-4") ----
solve_amg(SymmetricBlockGaussSeidelIterative(1), 2, "lin_SBGSIt1_2")
end

--solve_amg(ILUTScalar(1e-8), 2, "lin_ILUT_1e-8")


convCheck:set_maximum_steps(10)
convCheck:set_minimum_defect(1e-25)
convCheck:set_reduction(1e-26)




--solve_amg(CG(), GaussSeidel(), 2, "CG_GS_2")
util.PrintTable(results)
--util.printFormattedTable({results["lin_GS_2"],results["lin_GS_4"],results["lin_GS_10"]}, {"lastReduction", "time"})   


do return end

print("Apply solver: LINEAR SOLVER.")
	print("Init solver for operator.")
	linSolver = LinearSolver()
	if false then
		ilut = ILUTScalar(0.0)
		ilut:set_info(true)
		convCheck:set_maximum_steps(util.GetParamNumber("-sm"))
		linSolver:set_preconditioner(BlockGaussSeidel(1))
		--linSolver:set_preconditioner(GaussSeidel())
		--linSolver:set_preconditioner(Jacobi(0.6))
	else
		linSolver:set_preconditioner(amg)
	end
	
	linSolver:set_convergence_check(convCheck)
	linSolverA = linSolver
	linSolverA:init(A)
	print(amg:level_information_str())
	
	srand(iSeed)
	b:assign(b2)
	u:set_random(-1.0, 1.0)	
	domainDisc:adjust_solution(u)
	tBefore = os.clock()
	linSolver:apply_return_defect(u,b)
	tSolve = os.clock()-tBefore
	ScaleDomain(dom, MakeVec(0,0,0), MakeVec(1/stretchX, 1/stretchY, 1.0))
	SaveVectorForConnectionViewer(u, "amg_u_linear.vec")
	print("done, took "..tSolve.." s")

	lastReduction = convCheck:defect()/convCheck:previous_defect()
	avgReduction = convCheck:avg_rate()
	steps = convCheck:step()
	do return end
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
	print("done, took "..tSolveCG.." s")
	
	CGlastReduction = convCheck:defect()/convCheck:previous_defect()
	CGavgReduction = convCheck:avg_rate()
	CGsteps = convCheck:step()







-- ugshell -ex amg/famgLEcheck.lua -numRefs 5 -maxLevels 4 -stretchX 20  -delta 0.1 -Hinv 0 -maxParents 3  -threePenalty 1.4 --> 1.4



if ProcRank() == 0 then
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

