-------------------------------------------------------------------------------
--[[!
-- \file amg/laplace_sin_filtering.lua
-- \{
-- We set a Laplace-problem up here
--
]]--
--------------------------------------------------------------------------------

SetMinSecondsUntilProgress(3)
--PrintBuildConfiguration()

ug_load_script("setup.lua")
outfile = util.GetParam("-outfile", nil, "output directory")

if true then
	amgdisc.LoadStdDomAndApproxSpace()
else

	gridName = util.GetParam("-grid", "../skin/2D/geometry/brick/skin2d.aniso.ugx")
	requiredSubsets = {"LIP", "COR", "BOTTOM_SC", "TOP_SC"}
	dom = amgdisc.MyLoadGrid({ dim=2, gridName = gridName, requiredSubsets = requiredSubsets} )
	boundaries = {"TOP_SC", "BOTTOM_SC"}

	approxSpace = amgdisc.CreateStdApproxSpace(dom)
	domainDisc = amgdisc.CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries,
		diffusion={LIP=1, COR=10}
		})
end


--if NumProcs() <= 16 then
	--LogAllOutput("/Users/mrupp/")
--end
	
--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
precond = Jacobi()
OrderCuthillMcKee(approxSpace, true)

--------------------------------------------------------------------------------
-- Solve
--------------------------------------------------------------------------------

-- 1. Create Algebra
write("Assembling...")
A, u, b = amgdisc.AssembleAxb(approxSpace, domainDisc)


--	ScaleDomain(dom, MakeVec(0,0,0), MakeVec(1/stretchX, 1/stretchY, 1.0))

if numRefs >= 0 then
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(util.GetBestLU())
	gmgSmoother = SymmetricGaussSeidel()
	--Jacobi(0.5)
	
	gmg:set_smoother(gmgSmoother);
	
	gmg:set_cycle_type(1)
	gmg:set_rap(false)

	--gmg = Jacobi(0.66)
	gmgConvCheck = ConvCheck(30, 1e-12, 1e-12)
	linSolverGMG = LinearSolver(gmg, gmgConvCheck)
		
	print(linSolverGMG:config_string())
	
	--for v, i in pairs({1, 2, 4, 8}) do
	for v, i in pairs({1}) do
	gmg:set_num_presmooth(i)
	gmg:set_num_postsmooth(i)
	linSolverGMG:init(A)
	gmgRes=solve_random_d(linSolverGMG, approxSpace, gmgConvCheck)
	end
		
	print(gmgRes)
	
	if util.HasParamOption("-writeMatrices") then
		SaveVectorForConnectionViewer(gmgRes.u, "resGMG.vec")		
	end	
end


amg = FAMGPreconditioner()
amgsetup.SetFAMGCommandLineParameters(amg, {
	smooths=4,
	maxNodesForBase=2000,
	delta=0.1	
	})

if util.HasParamOption("-beam") then
print("BEAM")
	amg:set_beam(true)
	amg:set_external_coarsening(true)
	
end


sgs = SymmetricGaussSeidel()
amg:set_presmoother(sgs)
amg:set_postsmoother(sgs)

--SetDebugLevel(debugID.AMG.Coarsening, 4)


if util.HasParamOption("-writeMatrices") then
	matrixWritePath=CreateOutputDirectory("RAMmatrices")
	--matrixWritePath = "/Users/mrupp/matrices/"
	amg:set_matrix_write_path(matrixWritePath)
	amg:set_write_f_values(true)
	--amg:set_debug_write_all_vectors(true)
end
amg:set_write_testvectors(true)

-- 2: start with best two, increase : best opComplexity, also best convergence?
-- 3: start with best 6, decrease as long neg
-- 4: start with 8, reduce as long better
-- 3 doch besser?
amg:set_xc_calc_interpolation_routine(2)
amg:set_xc_calc_interpolation_routine(3)
amg:set_xc_calc_interpolation_routine(3)

-- unterschied sichtbar bei  ugshell -ex amg/famg_laplace.lua -33k -numRefs 0 -amgMaxLevels 5 -writeMatrices -amgTV lin -amgXC
--amg:set_xc_calc_interpolation_routine(4)

amg:set_deficit_solve_tolerance(1e-20)


--SetDebugLevel(debugID.AMG.Coarsening, 5)
--SetDebugLevel(debugID.AMG, 5)
amg:set_XC_LE_griebel(true)

--amg:add_watchnode(0, 1453)
--amg:add_watchnode(0, 175)


--amg:set_full_block_jacobi(true)
amg:set_fsmoothing(true)
amg:set_strong_connection_external(0.1)

baseEps = util.GetParamNumber("-baseEps", 0.0)
if baseEps > 0.0 then
	baseSolver = LinearSolver(Jacobi(0.8), ConvCheck(100000, 1e-20, baseEps, false, false))
	amg:set_base_solver(baseSolver)
end 

convCheck = ConvCheck(30, 1e-20, 1e-20)
--convCheck:set_linear_operator(A)

--SetDebugLevel(debugID.AMG, 5)
--SetDebugLevel(debugID.AMG.FAMG.Interpolation, 0)
amg:set_verbose(5)

	
--SetDebugLevel(debugID.AMG.Agglomeration, 5)
amg:set_min_nodes_on_one_processor(400)
amg:set_preferred_nodes_on_one_processor(600)
--amg:set_debug_write_all_vectors(true)

if util.HasParamOption("-posamg") then
	amg = POSAMGPreconditionerCPU1()
	amg:set_matrix_write_path(CreateOutputDirectory("RAMmatrices"))
	amg:set_position_provider( GridFunctionPositionProvider(GridFunction(approxSpace)) )
end
linSolver = LinearSolver(amg, convCheck)
solver = linSolver

print(solver:config_string())
print("Init Solver...")

solver:init(A)

lis=amg:level_information_str()
	print(lis)
print("Solve...")


--solve_random(solver, approxSpace)

amg:set_one_init(true)

if util.HasParamOption("-beam") then
	amg:set_nr_of_preiterations_at_check(20)
	amg:set_fsmoothing(false)
	amg:set_external_coarsening(true)
	sm = SparseBlockGaussSeidel2(3)
	bgsRes = test_smoothsteps(sm, convCheck, {1, 2, 4, 16, 64})
	for i, v in pairs(bgsRes) do
		SaveVectorForConnectionViewer(bgsRes[i].linV.u, "res"..i..".vec")
	end
elseif true then
	--sgsRes = test_smoothsteps(SymmetricGaussSeidel(), convCheck, {1, 2, 4, 8, 16, 32, 64, 128})
	--sgsRes = test_smoothsteps(SymmetricGaussSeidel(), convCheck, {1, 2, 4, 8, 16})
	--sgsRes = test_smoothsteps(SymmetricGaussSeidel(), convCheck, {1, 2, 4})
	sgsRes = test_smoothsteps(SymmetricGaussSeidel(), convCheck, {4})
	--EconvCheck = EnergyConvCheck(50, 1e-20, 1e-20)
	--EconvCheck:set_linear_operator(A)
	--sgsRes = test_smoothsteps(SymmetricGaussSeidel(), EconvCheck, {2})
	--jacRes = test_smoothsteps(Jacobi(0.5), convCheck, {1, 2, 4, 8, 16})
	
	--jacRes = test_smoothsteps(Jacobi(0.6), {1, 2, 4, 8})
	if util.HasParamOption("-writeMatrices") then
		for i, v in pairs(sgsRes) do
			if sgsRes[i].linV then SaveVectorForConnectionViewer(sgsRes[i].linV.u, "res"..i..".vec") end			
			if jacRes then SaveVectorForConnectionViewer(jacRes[i].linV.u, "jacres"..i..".vec") end
		end
	end
	--jacRes = test_smoothsteps(Jacobi(0.6), {1, 2, 4, 8})	
end

if true then
	convCheck = ConvCheck(100, 1e-12, 1e-14, false)
	print("calculating smooth vector...")
	if util.HasParamOption("-beam") then
		res = solve_lin(LinearSolver(), SparseBlockGaussSeidel2(2), convCheck, "bgs")
	else
		res = solve_lin(LinearSolver(), Jacobi(0.6), convCheck, "jac")
	end
	if util.HasParamOption("-writeMatrices") then
		SaveVectorForConnectionViewer(res.u, "smooth.vec")
		SaveVectorForConnectionViewer(res.d, "smooth_d.vec")
	end
	print("done\n")
end


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
	{"commandline", util.GetCommandLine() }
	}
	
	table.append(stats, amgsetup.getAMGStats(amg))
	 
	if gmgRes ~= nil then util.flattenTable(gmgRes, "gmg", stats) end
	if sgsRes ~= nil then util.flattenTable(sgsRes, "sgs", stats) end
	if jacRes ~= nil then util.flattenTable(jacRes, "jac", stats) end
	if bgsRes ~= nil then util.flattenTable(bgsRes, "bgs", stats) end	
	
	util.flattenTable(getAMGLevelTable(amg), "level", stats)
	
	util.printStats(stats)
	if outfile ~= nil  then	
		util.writeFileStats(stats, outfile)
	end

end
	
WriteProfileData("pd.pdxml")

amg:check_vector(TestvectorsConstant(), "constant")
if dim == 2 then amg:check_vector(TestvectorsLinear2d(), "tvLin")
else amg:check_vector(TestvectorsLinear3d(), "tvLin") end


if false and NumProcs() == 1 then
	amg:set_preiterations_mimum_defect_at_check(1e-30)	
	amg:set_nr_of_preiterations_at_check(20)
	amg:set_num_postsmooth(4)
	amg:set_num_presmooth(4)
	amg:set_postsmoother(SymmetricGaussSeidel())
	amg:set_presmoother(SymmetricGaussSeidel())
	amg:init_smoothers()
	
	--for i, v in pairs({1, 2, 4, 8, 16, 32, 64, 128}) do
	--for i, v in pairs({4}) do
	--	amg:set_num_postsmooth(i)
		--amg:set_num_presmooth(i)
		amg:check()
	--end
end
