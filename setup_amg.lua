amgsetup = {}



--! @return a table which represents a mapping name -> smoother
--! can be used in util.GetParamFromList 
function amgsetup.getSmootherList()
	return {
		jac05 = Jacobi(0.5),
		jac06 = Jacobi(0.6),
		jac07 = Jacobi(0.7),
		jac08 = Jacobi(0.8),
		jac09 = Jacobi(0.9),
		jac10 = Jacobi(1.0),
		ilu = ILU(),
		gs = GaussSeidel(),
		sgs = SymmetricGaussSeidel(),
		bgs = BackwardGaussSeidel()
	}	
end


--! @return a table which represents a mapping name -> ITestvectorCalculator 
--! can be used in util.GetParamFromList 
--! @sa amgsetup.GetFAMGTestvectorFromCommandLine
function amgsetup.GetFAMGTestvectorList(dim)
	if dim == 2 then
		return	
		{	const = TestvectorsConstant(), 
			lin = TestvectorsLinear2d(), 
			le = TestvectorsLinearElasticity2d() }		
	elseif dim == 3 then
		return
		{ 	const = TestvectorsConstant(),
			lin = TestvectorsLinear3d(),
			le = TestvectorsLinearElasticity3d() }
	else
		ug_assert(false, "dim = "..tostring(dim).." not supported")
	end
end

--! checks command line parameter -amgTV.
--! valid options are const, lin, le
--! default parameter is const
function amgsetup.GetFAMGTestvectorFromCommandLine(dim)
	return util.GetParamFromList("-amgTV", "const", amgsetup.GetFAMGTestvectorList(dim))
end

--! set FAMG amg parameters. if parameters are not given in p,
--! command line parameters or default values are chosen.
--! valid elements of p and their default values:
--! smooths=2 (-amgSmooths)
--! maxLevels=30 (-amgMaxLevels)
--! maxNodesForBase=1024 (-amgMaxNodesForBase)
--! bTwoLevel=false (-amgTwoLevel (option))
--! wcycle=false (-amgWcycle (option))
--! cycleType=1 (-amgCycle)
--! iTestvectorSmooths=0 (-amgTestvectorsSmooths)
--! testVectorSmoother=jac06 (-amgTVSmoother)
--! maxParents=2 (-amgMaxParents)
--! XC=false (-amgXC (option))
--! indirectInterpolation=false (-amgIndirectInterpolation)
function amgsetup.SetFAMGCommandLineParameters(amg, p)
	if false then amg = FAMGPreconditioner() end -- for UGIDE
	
	 p.smooths = util.GetParamNumber("-amgSmooths", p.smooths or 2)
	 p.maxLevels = util.GetParamNumber("-amgMaxLevels", p.maxLevels or 30)
	 p.maxNodesForBase = util.GetParamNumber("-amgMaxNodesForBase", p.maxNodesForBase or 1024)
	 p.bTwoLevel = p.bTwoLevel or util.HasParamOption("-amgTwoLevel")
	 if p.bTwoLevel then
	 	p.maxNodesForBase = 1
	 	p.maxLevels = 2
	 end
	
	if p.wcycle or util.HasParamOption("-amgWcycle") then
		if util.GetParamNumber("-amgCycle", "n/a") ~= "n/a" then 
			ug_warning("WARNING: Options -amgWcycle and -amgCycle affect the same thing! Only use one. WARNING")
		end  
		amg:set_cycle_type(2)
	else
		amg:set_cycle_type(1)
	end
	local amgCycleType = util.GetParamNumber("-amgCycle", p.cycleType or "n/a")
	if amgCycleType ~= "n/a" then
		amg:set_cycle_type(amgCycleType)
	end
		
	amg:set_testvectors(amgsetup.GetFAMGTestvectorFromCommandLine(dim))
	if util.HasParamOption("-constIfFailed") then
		amg:set_use_constant_if_failed(true)
	end
	
	
				
	p.smoother = util.GetParamFromList("-amgSmoother", "sgs", amgsetup.getSmootherList() or p.smoother)
	
	p.iTestvectorSmooths = util.GetParamNumber("-amgTestvectorSmooths", p.iTestvectorSmooths or 0)
	amg:set_testvector_smooths(p.iTestvectorSmooths)
	p.testVectorSmoother = util.GetParamFromList("-amgTVSmoother", "jac06", amgsetup.getSmootherList() or p.testVectorSmoother)
	amg:set_testvector_smoother(p.testVectorSmoother)
	amg:set_deficit_solve_tolerance(1e-12)

	amg:set_max_levels(p.maxLevels)
	amg:set_max_nodes_for_base(p.maxNodesForBase)
	amg:set_base_solver(util.GetBestLU())
		
		
	amg:set_num_postsmooth(p.smooths)
	amg:set_num_presmooth(p.smooths)	
	amg:set_max_parents(util.GetParamNumber("-amgMaxParents", p.maxParents or 2))
	
	amg:set_delta(util.GetParamNumber("-amgDelta", p.delta or 0.1))
	
	if p.XC or util.HasParamOption("-amgXC") then amg:set_external_coarsening(true) end
	if p.N2Interpolation or util.HasParamOption("-amgN2Interpolation") then amg:set_N2_interpolation(true) end
	
	
	uuu = GridFunction(approxSpace)
	uuuu = GridFunctionPositionProvider(uuu)
	amg:set_position_provider( uuuu )
	amg:set_indirect_interpolation(p.indirectInterpolation or util.HasParamOption("-amgIndirectInterpolation")) 	
 end
 
 function getAMGLevelTable(amg)
	local res = {}
	local levs = amg:get_used_levels()-1
	if levs > 7 then
		levs = 7
	end
	for i=0,levs do
		li = amg:get_level_information(i)
		res[i+1] = li:get_nr_of_nodes()
	end
	for i=levs+1, 7 do
		res[i+1] = 0
	end
	return res

end
 
 
function amgsetup.getAMGStats(amg)
	return 
	{ { "ndofs", amg:get_level_information(0):get_nr_of_nodes() },
	{ "tSetupAmg [s]", amg:get_timing_whole_setup_ms()/1000},
	
	{ "amg_base_preprocess [ms]", GetProfileNode("amg_base_preprocess"):get_avg_total_time_ms()},
	{ "amg_createDirectSolver [ms]", GetProfileNode("amg_createDirectSolver"):get_avg_total_time_ms()},
		
	{ "c_A", amg:get_operator_complexity()},
	{ "c_G", amg:get_grid_complexity()},
	{ "used Levels", amg:get_used_levels()},
	{ "maxBase", amg:get_max_nodes_for_base() },
	{ "delta", amg:get_delta() },
	{ "tv", util.GetParam("-amgTV", "const") }
	}	
end
 function amgsetup.getAMGStats(amg) return {} end