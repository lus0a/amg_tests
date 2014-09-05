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


amgsetup.LoadStdDomAndApproxSpace()

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")
--------------------------------------------------------------------------------
-- Solve
--------------------------------------------------------------------------------

-- 1. Create Algebra
A, u, b = AssembleAxb(approxSpace, domainDisc)

if true then
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(GetBestLU())
	
	gmgSmoother = util.GetParamFromList("-smoother", "sgs", amgsetup.getSmootherList())
	numSmooths = util.GetParamNumber("-numSmooths", 2)
		
	gmg:set_smoother(gmgSmoother);
	gmg:set_num_presmooth(numSmooths)
	gmg:set_num_postsmooth(numSmooths)
	gmg:set_cycle_type(1)
	
	print(gmg:config_string())

	--gmg = Jacobi(0.66)
	gmgConvCheck = ConvCheck(30, 1e-12, 1e-12)
	linSolverGMG = LinearSolver(gmg, gmgConvCheck)
	linSolverGMG:init(A)
	gmgRes=solve_random_d(linSolverGMG, approxSpace, gmgConvCheck)
	print(gmgRes)
	
	if util.HasParamOption("-writeMatrices") then
		SaveVectorForConnectionViewer(gmgRes.u, "resGMG.vec")
	end	
end
	


if ProcRank() == 0 then
	stats = {
	{ "order", util.GetParamOption("-order", 0)},	
	{ "procs", NumProcs() },
	--{ "numPreRefs", numPreRefs},
	{ "numRefs", numRefs },
	{ "nu", util.GetParamNumber("-nu")},
	{ "numSmooths", numSmooths},
	{ "dim", dim},
	{ "gridName", gridName},
	{ "main [ms]", GetProfileNode("main"):get_avg_total_time_ms()},
	{ "perform_refinement [ms]", GetProfileNode("perform_refinement"):get_avg_total_time_ms()},
	 
	{ "date", os.date("y%Ym%md%d") },
	{ "SVN Revision", GetSVNRevision()},
	{"host", GetBuildHostname()},
	{"commandline", util.GetCommandLine() }
	}
	
	if gmgRes ~= nil then flattenTable(gmgRes, "gmg", stats) end
	
	util.printStats(stats)
	if outfile ~= nil  then	
		util.writeFileStats(stats, outfile)
	end

end
	
WriteProfileData("pd.pdxml")
