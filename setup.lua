ug_load_script("ug_util.lua")
--ug_load_script("solver_util/setup_rsamg.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("util/mac_util.lua")


function SetDefaultParameters(tableParameter, default)
	for i, v in pairs(default) do
		if tableParameter[i] == nil then tableParameter[i] = v end
	end
end

function NeededParameters(p, parameters)
	for i, v in ipairs(parameters) do
		ug_assert(p.v ~=nil)	
	end
end



function SetFAMGCommandLineParameters(amg, p)
	if false then amg = FAMGPreconditioner() end
	 p.smooths = util.GetParamNumber("-amgSmooths", p.smooths or 2)
	 p.maxLevels = util.GetParamNumber("-amgMaxLevels", p.maxLevels or 8)
	 p.maxNodesForBase = util.GetParamNumber("-amgMaxNodesForBase", p.maxNodesForBase or 1024)
	 p.bTwoLevel = p.bTwoLevel or util.HasParamOption("-amgTwoLevel")
	 if p.bTwoLevel then
	 	p.maxNodesForBase = 1
	 	p.maxLevels = 2
	 end
	
	if p.wcycle or util.HasParamOption("-amgWcycle") then
		amg:set_cycle_type(2)
	else
		amg:set_cycle_type(1)
	end
		 
	p.iTestvectorSmooths = util.GetParamNumber("-amgTestvectorSmooths", p.iTestvectorSmooths or 4)
	amg:set_testvector_smooths(p.iTestvectorSmooths)
	amg:set_deficit_solve_tolerance(1e-12)

	amg:set_max_levels(p.maxLevels)
	amg:set_max_nodes_for_base(p.maxNodesForBase)
	if HasClassGroup("SuperLU") then amg:set_base_solver(SuperLU()) end	
		
	amg:set_num_postsmooth(p.smooths)
	amg:set_num_presmooth(p.smooths)
	amg:set_max_parents(util.GetParamNumber("-amgMaxParents", p.maxParents or 2))
	
	if p.XC or util.HasParamOption("-amgXC") then amg:set_external_coarsening(true) end
	
	uuu = GridFunction(approxSpace)
	uuuu = GridFunctionPositionProvider(uuu)
	amg:set_position_provider( uuuu )
	amg:set_indirect_interpolation(false) 	
 end

function last_random_conv(solver, approxSpace)
	local f = GridFunction(approxSpace)
	f:set_random(-1, 1)		
	solver:apply_return_defect(f, ZeroFunction())
	local lastReduction= convCheck:defect()/convCheck:previous_defect()
	print("solve random lastReduction = ".. lastReduction )
	return lastReduction
end

function solve_random(solver, approxSpace)
	local avgRate
	local lastReduction
	local convMsg
  	local res = solve_random_d(solver, approxSpace)
  	if res.avgRate < 1e-12 then
		convMsg = res.avgRate.." - EXACT!"
	else
		convMsg = res.avgRate
	end	
  	
	print("solve random avgConv = "..convMsg.." lastReduction = ".. res.lastReduction )
	if false and res.avgRate < 1e-12 then return "EXACT", "EXACT"
	else 
	return string.format("%2.2e", res.avgRate), string.format("%2.2e", res.lastReduction) 
	end
end

function solve_random_d(solver, approxSpace)
	local tBeginS = GetClockS()
	local f = GridFunction(approxSpace)
	f:set_random(-1, 1)		
	solver:apply_return_defect(f, ZeroFunction())
	local avgRate = convCheck:avg_rate()
	local lastReduction = convCheck:defect()/convCheck:previous_defect()
	local timeSolveS = GetClockS()-tBeginS
	return {avgRate=avgRate, lastReduction=lastReduction, timeSolveS=timeSolveS}	
end

function MyLoadGrid(p)
	SetDefaultParameters(p, 
		{dim=2, gridName="unit_square_01/unit_square_01_quads_2x2.ugx",
		requiredSubsets = {"Inner", "Boundary"}, cpu=1, dim=2
		} )
	
	-- refinements:
	numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
	numRefs    = util.GetParamNumber("-numRefs",    3, "number of refinements")
	
	if numPreRefs > numRefs then numPreRefs = numRefs end
		
	-- Display parameters (or defaults):
	print("CommandLine: " .. util.GetCommandLine())
	print("Parameters: ")
	s = StringTable()
	s:set(0, 0, "grid = " .. gridName)
	s:set(0, 1, "numRefs = " .. numRefs)
	s:set(1, 0, "numPreRefs = " .. numPreRefs)

	
	-- choose algebra
	InitUG(p.dim, AlgebraType("CPU", p.cpu));	
	
	if withRedist == true then
		dom = util.CreateDomain(gridName, 0, p.requiredSubsets)
		balancer.RefineAndRebalanceDomain(dom, numRefs)
	else
		dom = util.CreateAndDistributeDomain(p.gridName, numRefs, numPreRefs, p.requiredSubsets)
	end	
		
	stretchX = util.GetParamNumber("-stretchX", 1.0)
	stretchY = util.GetParamNumber("-stretchY", 1.0)
	print("scaling domain by ("..stretchX..", "..stretchY..")") 
	ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))
	
	distortDomain = util.GetParamNumber("-distortDomain", 0.0)
	if distortDomain > 0.0 then
		distortDomain = math.pow(d, numRefs+1)
		print("distortDomain ="..distortDomain)
		RandomizeDomain(dom, distortDomain, distortDomain, distortDomain)
	end			
		
	return dom
end

function CreateOutputDirectory(name, mbsize)
	local thePath
	name = "/"..name.."/"
	if GetOperatingSystem() == "apple" then
		if mbsize == nil then mbsize = 1024 end
		thePath = util.createMacRAMDisk("ug4RAMDisk", mbsize)..name		
	else
		thePath = os.getenv("HOME")..name
	end
	CreateDirectory(thePath) 
	return thePath		
end

function LogAllOutput(thePath)
	a = GetLogAssistant()
	a:set_output_process(-1)
	if paramOutputProcess ~= -1 and ProcRank() ~= paramOutputProcess then a:enable_terminal_output(false) end
	a:enable_file_output(true, thePath.."log_p"..ProcRank()..".txt")	
end


function CreateStdApproxSpace(dom)
	-- create Approximation Space
	approxSpace = ApproximationSpace(dom)
	approxSpace:add_fct("c", "Lagrange", 1)
	
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	-- lets order indices using Cuthill-McKee
	--OrderCuthillMcKee(approxSpace, true);
	return approxSpace
end

function CreateStdDisc(p)
	local dim = approxSpace:get_dim()
	 
	local elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
	elemDisc:set_diffusion(1.0)
	--elemDisc:set_diffusion("mydiffusion")
	elemDisc:set_source("zero_"..dim)
		
	local dirichletBND = DirichletBoundary()
	for i,v in pairs(p.boundaries) do
		dirichletBND:add("dirichlet_zero_"..dim, "c", v)
	end
	
	local domainDisc = DomainDiscretization(p.approxSpace)
	domainDisc:add(elemDisc)
	domainDisc:add(dirichletBND)
	
	return domainDisc
end

function CreateLEDisc(parameters)
	local dom = parameters.Domain
	local BoundaryName = parameters.BoundaryName
	local dim = dom:get_dim()
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
	
	approxSpace = ApproximationSpace(dom)
		
	-- Create, Load, Refine and Distribute Domain
		
	-----------------------------------------------------------------
	--  Setup FV Convection-Diffusion Element Discretization
	-----------------------------------------------------------------
	
	if bLinearElasticity then
		approxSpace:add_fct("ux", "Lagrange", 1)
		approxSpace:add_fct("uy", "Lagrange", 1)
		approxSpace:init_levels()
		approxSpace:init_top_surface()
		approxSpace:print_statistic()
		
		if dim == 2 then
			elemDisc = SmallStrainMechanics("ux, uy", "Inner")
		else
			elemDisc = SmallStrainMechanics("ux, uy, uz", "Inner")
		end
		matLaw = HookeLaw()
		nu = util.GetParamNumber("-nu", 0.4)
		matLaw:set_hooke_elasticity_tensor_E_nu(20.6316, nu)
		elemDisc:set_material_law(matLaw)
	else
		approxSpace:add_fct("c", "Lagrange", 1)          
		approxSpace:init_levels()
		approxSpace:init_top_surface()
		elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
		elemDisc:set_diffusion(1)
		elemDisc:set_source(0)	
	end
	-----------------------------------------------------------------
	--  Setup Dirichlet Boundary
	-----------------------------------------------------------------

	dirichletBND = DirichletBoundary()
	function setBoundary(boundaryName)
		if bLinearElasticity then
			dirichletBND:add("dirichlet_zero_"..dim, "ux", boundaryName)
			dirichletBND:add("dirichlet_zero_"..dim, "uy", boundaryName)
			if dim == 3 then
				dirichletBND:add("dirichlet_zero_"..dim, "uz", boundaryName)
			end
		else
			dirichletBND:add("dirichlet_zero_"..dim, "c", boundaryName)
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
			end
		end
	end
	-------------------------------------------
	--  Setup Domain Discretization
	-------------------------------------------
	
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(elemDisc)
	domainDisc:add(dirichletBND)
	-- domainDisc:add(neumannDisc)
	
	-- ? not necessary
	--OrderCuthillMcKee(approxSpace, false);
	return approxSpace, domainDisc
end

function zero_2(x, y, t)
	return 0;
end

function zero_3(x, y, z, t)
	return 0;
end

function dirichlet_zero_2(x, y, t)
	if x == 0 then
		return false, 1
	else
		return true, 0
	end
end

function dirichlet_zero_3(x, y, z, t)
	return true, 0
end

function dirichlet_one_2(x, y, t)
	return true, 1
end

function dirichlet_one_3(x, y, z, t)
	return true, 1
end

function mydiffusion(x, y, t)
	return 1,0,0,0.01
end


function smooth(gf, smoother, n, A)
	local convCheck = ConvCheck(n, 0, 0, false, true)
	local a = LinearSolver(smoother, convCheck)
	a:init(A)
	a:apply(gf, ZeroFunction())
	return convCheck:defect()/convCheck:previous_defect()
end


function CreateGridFunction(func, fct)
	if fct == nil then fct = 0 end
	local gf = GridFunction(approxSpace)
	Interpolate(func, gf, approxSpace:name(fct))
	--smooth(gf, Jacobi(0.6), 4)
	return gf
end



function ZeroFunction()
	local gf = GridFunction(approxSpace)
	gf:set(0.0)
	return gf
end

function RandomFunction()
	local gf = GridFunction(approxSpace)
	gf:set_random(-1.0, 1.0)
	return gf
end


function CreateSinGridFunction(x, y)
	nuX = x
	nuY = y	
	local a= CreateGridFunction("myEVFunction")
	domainDisc:adjust_solution(a)
	--VecScaleAssign(a, 1/ VecNorm(a), a)
	return a
end


function CreateSinGridFunction2(x1, y1, x2, y2)
	local gf = GridFunction(approxSpace)
	nuX = x1
	nuY = y1
	SetGridFunction(gf, "myEVFunction", 0)
	nuX = x2
	nuY = y2
	SetGridFunction(gf, "myEVFunction", 1)
	domainDisc:adjust_solution(gf)
	return gf
end


function SetGridFunction(gf, func, fct)	
	Interpolate(func, gf, approxSpace:name(fct))	
	return gf
end

function myEVFunction2(x, y, t)
	return math.sin(math.pi*nuX*x)*math.sin(math.pi*y*nuY)
end



function myEVFunction(x, y, t)
	if nuX==0 or nuY==0 then return (1+math.sin(math.pi*nuX*x))*(1+math.sin(math.pi*y*nuY)) end
	return math.sin(math.pi*nuX*x)*math.sin(math.pi*y*nuY)
end


function flattenTable(tab, name, res)
	
	local tName = ""
	if name ~= nil then
		tName = name.."."
	end
	if type(tab) == "table" then
		for i, v in pairs(tab) do
			flattenTable(v, tName..i, res)
		end
		return res
	else
		table.insert(res, {name, tab})
	end
end
	
function getAMGStats(amg)
	return 
	{ { "ndofs", amg:get_level_information(0):get_nr_of_nodes() },
	{ "tSetupAmg [s]", amg:get_timing_whole_setup_ms()/1000},
	
	{ "amg_base_preprocess [ms]", GetProfileNode("amg_base_preprocess"):get_avg_total_time_ms()},
	{ "amg_createDirectSolver [ms]", GetProfileNode("amg_createDirectSolver"):get_avg_total_time_ms()},
		
	{ "c_A", amg:get_operator_complexity()},
	{ "c_G", amg:get_grid_complexity()},
	{ "used Levels", amg:get_used_levels()},
	{ "maxBase", amg:get_max_nodes_for_base() },
	{ "delta", amg:get_delta() }
	}	
end

function solve_lin(linSolver, amg, name)
	linSolver:set_convergence_check(convCheck)	
	linSolver:set_preconditioner(amg)
	linSolver:init(A)
	return solve_random_d(linSolver, approxSpace)
	
end

function solve_amg(smoother, numSmooths, name)


	amg:set_presmoother(smoother)
	amg:set_postsmoother(smoother)
	amg:set_num_presmooth(numSmooths)
	amg:set_num_postsmooth(numSmooths)
	
	cgSolver = BiCGStab()
	cgSolver:set_preconditioner(amg)	
	amg:init_smoothers()
	
	print(name..": "..numSmooths.."x "..smoother:config_string())
		
	print("CG, "..name)
	local ret = {}
	ret.cg = solve_lin(CG(), amg, "CG("..name..")")
	print("LinearSolver, "..name)
	ret.lin = solve_lin(LinearSolver(), amg, "LS("..name..")")
	return ret
end
