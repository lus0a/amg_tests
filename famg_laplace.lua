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

InitSignals()

outfile = util.GetParam("-outfile", nil, "output directory")


if true then
	amgsetup.LoadStdDomAndApproxSpace()
else

	gridName = util.GetParam("-grid", "../skin/2D/geometry/brick/skin2d.aniso.ugx")
	requiredSubsets = {"LIP", "COR", "BOTTOM_SC", "TOP_SC"}
	dom = MyLoadGrid({ dim=2, gridName = gridName, requiredSubsets = requiredSubsets} )
	boundaries = {"TOP_SC", "BOTTOM_SC"}

	approxSpace = CreateStdApproxSpace(dom)
	domainDisc = CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries,
		diffusion={LIP=1, COR=10}
		})
end
OrderCuthillMcKee(approxSpace, true)
	

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

convCheck = ConvCheck(5, 1e-5, 1e-5, false)
precond = Jacobi()


--------------------------------------------------------------------------------
-- Solve
--------------------------------------------------------------------------------

-- 1. Create Algebra
A, u, b = AssembleAxb(approxSpace, domainDisc)write("Assembling...")

amg = FAMGPreconditioner()
SetFAMGCommandLineParameters(amg, {
	--delta=0.030473,
	delta=0.1,
	smooths=4,
	maxNodesForBase=2000	
	})


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

amg:set_delta(0.03)
-- 2: start with best too, increase : best opComplexity, also best convergence?
-- 3: start with best 6, decrease as long neg
-- 4: start with 8, reduce as long better
amg:set_xc_calc_interpolation_routine(2)
amg:set_xc_calc_interpolation_routine(3)
-- unterschied sichtbar bei  ugshell -ex amg/famg_laplace.lua -33k -numRefs 0 -amgMaxLevels 5 -writeMatrices -amgTV lin -amgXC
--amg:set_xc_calc_interpolation_routine(4)

--amg:set_delta(0.045)
amg:set_delta(0.2)
amg:set_delta(util.GetParamNumber("-delta", amg:get_delta()))
amg:set_deficit_solve_tolerance(1e-15)

amg:set_XC_LE_griebel(true)



amg:set_fsmoothing(true)
amg:set_delta(0.001, 4)
--amg:set_verbose(5)
--amg:add_watchnode(2, 939)
--amg:add_watchnode(1, 2800)
--amg:add_watchnode(0, 1084)
-- amg:add_watchnode(1, 839)
--amg = RSAMGPreconditioner()

--amg:add_watchnode(1, 21066)

--amg = RSAMGPreconditioner()
amg:set_strong_connection_external(0.1)
convCheck = ConvCheck(30, 1e-12, 1e-14)
linSolver = LinearSolver(amg, convCheck)
solver = linSolver

print(solver:config_string())
print("Init Solver...")

--SetDebugLevel(debugID.AMG, 5)
--SetDebugLevel(debugID.AMG.FAMG.Interpolation, 0)
amg:set_verbose(5)


if util.HasParamOption("-beam") then
	amg:set_beam(true)
end
ScaleDomain(dom, MakeVec(0,0,0), MakeVec(1/stretchX, 1/stretchY, 1.0))




if true then
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(numRefs-1)
	gmg:set_base_solver(SuperLU())
	iluts = GaussSeidel()
		
	gmg:set_smoother(iluts);
	gmg:set_num_presmooth(32)
	gmg:set_num_postsmooth(32)
	gmg:set_cycle_type(2)

	
	linSolverGMG = LinearSolver(gmg, ConvCheck(30, 1e-12, 1e-12))
	linSolverGMG:init(A)
	resGMG=solve_random_d(linSolverGMG, approxSpace)
	print(resGMG)
	
	SaveVectorForConnectionViewer(resGMG.u, "resGMG.vec")
	do return end
end
	


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
	sm = SparseBlockGaussSeidel2(2)
	bgsRes = test_smoothsteps(sm, {1, 2, 4, 16, 64})
else
	sgsRes = test_smoothsteps(SymmetricGaussSeidel(), {1, 2, 4, 8, 16, 32, 64})
	sgsRes = test_smoothsteps(Jacobi(0.6), {1, 2, 4, 8, 16})
	
	--jacRes = test_smoothsteps(Jacobi(0.6), {1, 2, 4, 8})
	if false then
	for i, v in pairs(sgsRes) do
		SaveVectorForConnectionViewer(sgsRes[i].lin.u, "res"..i..".vec")
		if jacRes then SaveVectorForConnectionViewer(jacRes[i].lin.u, "jacres"..i..".vec") end
	end
	end
	--jacRes = test_smoothsteps(Jacobi(0.6), {1, 2, 4, 8})	
end

if true then
	convCheck = ConvCheck(100, 1e-12, 1e-14, false)
	print("calculating smooth vector...")
	if util.HasParamOption("-beam") then
		res = solve_lin(LinearSolver(), SparseBlockGaussSeidel2(2), "bgs")
	else
		res = solve_lin(LinearSolver(), Jacobi(0.6), "jac")
	end
	SaveVectorForConnectionViewer(res.u, "smooth.vec")
	SaveVectorForConnectionViewer(res.d, "smooth_d.vec")
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
	
	table.append(stats, getAMGStats(amg)) 
	
	if sgsRes ~= nil then flattenTable(sgsRes, "sgs", stats) end
	if jacRes ~= nil then flattenTable(jacRes, "jac", stats) end
	if bgsRes ~= nil then flattenTable(bgsRes, "bgs", stats) end	
	
	flattenTable(getAMGLevelTable(amg), "level", stats)
	
	util.printStats(stats)
	if outfile ~= nil  then	
		util.writeFileStats(stats, outfile)
	end

end
	
WriteProfileData("pd.pdxml")

amg:check_vector(TestvectorsConstant(), "constant")
if dim == 2 then amg:check_vector(TestvectorsLinear2d(), "tvLin")
else amg:check_vector(TestvectorsLinear3d(), "tvLin") end


if NumProcs() == 1 and true then
	amg:set_preiterations_mimum_defect_at_check(1e-9)
	amg:set_nr_of_preiterations_at_check(5)
	amg:set_num_postsmooth(4)
	amg:set_num_presmooth(4)
	--amg:set_postsmoother(SymmetricGaussSeidel())
	--amg:set_presmoother(SymmetricGaussSeidel())
	--amg:init_smoothers()
	
	--for i, v in pairs({1, 2, 4, 8, 16, 32, 64, 128}) do
	for i, v in pairs({4}) do
		amg:set_num_postsmooth(i)
		amg:set_num_presmooth(i)
		amg:check()
	end
end
