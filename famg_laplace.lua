-------------------------------------------------------------------------------
--[[!
-- \file amg/laplace_sin_filtering.lua
-- \{
-- We set a Laplace-problem up here
--
]]--
--------------------------------------------------------------------------------


--PrintBuildConfiguration()

ug_load_script("setup.lua")

InitSignals()

dim = util.GetParamNumber("-dim", 2)
outdir = util.GetParam("-outdir", nil, "output directory")

if dim == 3 then
gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_4x4x4.ugx")
dom = MyLoadGrid({ dim=3, gridName = gridName, requiredSubsets = requiredSubsets} )
else
gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")
dom = MyLoadGrid({ dim=2, gridName = gridName, requiredSubsets = requiredSubsets} )

end


requiredSubsets = {"Inner", "Boundary"}
boundaries = {"Boundary"}

util.CheckAndPrintHelp("");


approxSpace = CreateStdApproxSpace(dom)
domainDisc = CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries})

	

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

convCheck = ConvCheck(30, 1e-11, 1e-12, false)
precond = Jacobi()


--------------------------------------------------------------------------------
-- Solve
--------------------------------------------------------------------------------

-- 1. Create Algebra
A = MatrixOperator()
b = GridFunction(approxSpace)
u = GridFunction(approxSpace)


u:set(0.0)
domainDisc:adjust_solution(u)

-- 3. assemble matrix and rhs
write("Assembling...")
domainDisc:assemble_linear(A, b)
--SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
print(" done.")

amg = FAMGPreconditioner()
SetFAMGCommandLineParameters(amg, {
	--delta=0.030473,
	delta=0.1,
	smooths=4,
	maxNodesForBase=2000	
	})
amg:set_deficit_solve_tolerance(1e-15)

sgs = SymmetricGaussSeidel()
amg:set_presmoother(sgs)
amg:set_postsmoother(sgs)

matrixWritePath=CreateOutputDirectory("RAMmatrices")
amg:set_matrix_write_path(matrixWritePath)
amg:set_write_testvectors(true)

amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
amg:set_testvector_smooths(0)
amg:set_testvector_smoother(Jacobi(0.6))

amg:set_delta(0.1)
--amg:set_use_constant_if_failed(true)
amg:set_XC_LE_griebel(true)
amg:set_write_f_values(true)
--amg:set_debug_write_all_vectors(true)
amg:set_fsmoothing(true)
amg:set_testvectors(TestvectorsLinear2d())
--amg:set_testvectors(TestvectorsConstant())
--amg:set_verbose(5)
--amg:add_watchnode(0,0)
--amg:add_watchnode(0, 1084)
--amg:add_watchnode(1, 1084)
--amg = RSAMGPreconditioner()

--amg = RSAMGPreconditioner()
convCheck = ConvCheck(30, 1e-12, 1e-14)
linSolver = LinearSolver(amg, convCheck)
solver = linSolver

print(solver:config_string())
print("Init Solver...")

--SetDebugLevel(debugID.AMG, 5)
--SetDebugLevel(debugID.AMG.FAMG.Interpolation, 0)
amg:set_verbose(5)


solver:init(A)

lis=amg:level_information_str()
	print(lis)
print("Solve...")


--solve_random(solver, approxSpace)

amg:set_one_init(true)

sgs02 = solve_amg(SymmetricGaussSeidel(), 2, "SGSx02")
sgs32 = solve_amg(SymmetricGaussSeidel(), 32, "SGSx32")


if ProcRank() == 0 then
	stats = {
	
	{ "procs", NumProcs() },
	--{ "numPreRefs", numPreRefs},
	{ "numRefs", numRefs },
	{ "dim", dim},
	{ "gridName", gridName},
	{ "main [ms]", GetProfileNode("main"):get_avg_total_time_ms()},
	{ "perform_refinement [ms]", GetProfileNode("perform_refinement"):get_avg_total_time_ms()},
	 
	{ "date", os.date("y%Ym%md%d") },
	{ "SVN Revision", GetSVNRevision()},
	{"host", GetBuildHostname()},
	{"commandline", util.GetCommandLine() } }
	
	table.append(stats, getAMGStats(amg)) 
	
	flattenTable(sgs02, "sgs02", stats)
	flattenTable(sgs32, "sgs32", stats)
	
	flattenTable(getAMGLevelTable(amg), "level", stats)
	
	util.printStats(stats)
	if outdir ~= nil  then	
		util.writeFileStats(stats, outdir.."/stats.txt")
	end

end
	

