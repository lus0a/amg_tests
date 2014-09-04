ug_load_script("ug_util.lua")
--ug_load_script("solver_util/setup_rsamg.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("util/mac_util.lua")

EnableLUA2VM(true)
SetDebugLevel(debugID.LUACompiler, 4)
amgsetup = {}

function amgsetup.getSmootherList()
	return {
		sgs = SymmetricGaussSeidel(),
		jac05 = Jacobi(0.5),
		jac06 = Jacobi(0.6),
		jac07 = Jacobi(0.7),
		jac08 = Jacobi(0.8),
		jac09 = Jacobi(0.9),
		jac10 = Jacobi(1.0),
		ilu = ILU(),
		gs = GaussSeidel(),
		bgs = BackwardGaussSeidel() }	
end

function SetDefaultParameters(tableParameter, default)
	for i, v in pairs(default) do
		if tableParameter[i] == nil then tableParameter[i] = v end
	end
end

function NeededParameters(p, parameters)
	for i, v in pairs(parameters) do
		if p[v] == nil then
			print("\ntable ")
			print(p)
			print("needs parameter \""..v.."\" !")
			ug_assert(false, "table needs parameter \""..v.."\" !")
		end	
	end
end


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

function amgsetup.GetFAMGTestvectorFromCommandLine(dim)
	return util.GetParamFromList("-amgTV", "const", amgsetup.GetFAMGTestvectorList(dim))
end

function SetFAMGCommandLineParameters(amg, p)
	if false then amg = FAMGPreconditioner() end
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
	local amgCycleType = util.GetParamNumber("-amgCycle", "n/a")
	if amgCycleType ~= "n/a" then
		amg:set_cycle_type(amgCycleType)
	end
		
	amg:set_testvectors(amgsetup.GetFAMGTestvectorFromCommandLine(dim))
	if util.HasParamOption("-constIfFailed") then
		amg:set_use_constant_if_failed(true)
	end
				
	
	p.iTestvectorSmooths = util.GetParamNumber("-amgTestvectorSmooths", p.iTestvectorSmooths or 0)
	amg:set_testvector_smooths(p.iTestvectorSmooths)
	p.testVectorSmoother = util.GetParamFromList("-amgTVSmoother", "jac06", amgsetup.getSmootherList() or p.testVectorSmoother)
	amg:set_testvector_smoother(p.testVectorSmoother)
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
	amg:set_indirect_interpolation(util.HasParamOption("-amgIndirectInterpolation")) 	
 end

function last_random_conv(solver, approxSpace)
	return solve_random_d(solver, approxSpace).lastReduction
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

function GetConvCheckTable(convCheck, tBeginS)
	local avgRate = convCheck:avg_rate()
	local lastReduction = convCheck:defect()/convCheck:previous_defect()
	local timeSolveS
	if tBeginS then
		timeSolveS = GetClockS()-tBeginS
	else
		timeSolveS =nil
	end	
	return {avgRate=avgRate, lastReduction=lastReduction, timeSolveS=timeSolveS, steps=convCheck:step()}
end

function solve_random_d(solver, approxSpace)
	local tBeginS = GetClockS()
	local f = GridFunction(approxSpace)
	f:set_random(-1, 1)		
	local d = ZeroFunction()
	solver:apply_return_defect(f, d)
	local res= GetConvCheckTable(convCheck, tBeginS)
	res.u = f
	res.d = d
	return res		
end

function amgsetup.MyLoadGrid(p)
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
	write("Parameters: ")
	s = StringTable()
	s:set(0, 0, "grid = " .. gridName)
	s:set(0, 1, "numPreRefs = " .. numPreRefs)
	s:set(0, 2, "numRefs = " .. numRefs)
	print(s)

	
	-- choose algebra
	--InitUG(p.dim, AlgebraType("CPU", p.cpu));
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
		distortDomain = math.pow(distortDomain, numRefs+1)
		print("distortDomain ="..distortDomain)
		RandomizeDomain(dom, distortDomain, distortDomain, distortDomain)
	end			
		
	return dom
end

function getHome()
	local p=os.getenv("HOME")
	if p == nil then
		return "."
	else
		return p
	end
end

function CreateOutputDirectory(name, mbsize)
	local thePath
	name = "/"..name.."/"
	if GetOperatingSystem() == "apple" then
		if mbsize == nil then mbsize = 2048 end
		thePath = util.createMacRAMDisk("ug4RAMDisk", mbsize)..name		
	else
		thePath = getHome()..name
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
	NeededParameters(p, {"boundaries", "approxSpace"})

	local domainDisc = DomainDiscretization(p.approxSpace)
	local dim = p.approxSpace:get_dim()
	 
	 if p.diffusion == nil then
	 	local elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
		elemDisc:set_diffusion(1.0)
		--elemDisc:set_diffusion("mydiffusion")
		elemDisc:set_source("zero_"..dim)
		domainDisc:add(elemDisc)
	else
		for i, v in pairs(p.diffusion) do
			local elemDisc = ConvectionDiffusion("c", i, "fv1")
			elemDisc:set_diffusion(v)
			elemDisc:set_source("zero_"..dim)
			domainDisc:add(elemDisc)
		end
	end
		
	local dirichletBND = DirichletBoundary()
	
	for i,v in pairs(p.boundaries) do
		dirichletBND:add("dirichlet_zero_"..dim, "c", v)
	end
		
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
	
	InitUG(dim, AlgebraType("CPU", cpu))
	
	
	
	approxSpace = ApproximationSpace(dom)
		
	-- Create, Load, Refine and Distribute Domain
		
	-----------------------------------------------------------------
	--  Setup FV Convection-Diffusion Element Discretization
	-----------------------------------------------------------------
	
	approxSpace:add_fct("ux", "Lagrange", 1)
	approxSpace:add_fct("uy", "Lagrange", 1)
	if dim == 3 then
		approxSpace:add_fct("uz", "Lagrange", 1)
	end
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()
	
	
	varNames = "ux, uy"
	if dim == 3 then varNames = "ux, uy, uz" end	
	
		
	domainDisc = DomainDiscretization(approxSpace)
	if true then
		if parameters.le == nil then
			parameters.le = {}
			parameters.le["Inner"] = {}
			parameters.le["Inner"].nu = util.GetParamNumber("-nu", 0.4)
			parameters.le["Inner"].E = util.GetParamNumber("-leE", 20.6316)		
		end
		
		for i, v in pairs(parameters.le) do
			elemDisc = SmallStrainMechanics(varNames, i)	
			matLaw = HookeLaw()		
			matLaw:set_hooke_elasticity_tensor_E_nu(v.E, v.nu)
			elemDisc:set_material_law(matLaw)
			elemDisc:set_quad_order(2) --- !!!!!!!!
			print("--------------------------------------")
			print("In Subdomain "..i..":")
			print(elemDisc:config_string())
			domainDisc:add(elemDisc)
		end
		
	end
	if false then
		elemDisc = SmallStrainMechanics(varNames, "Inner")	
		matLaw = HookeLaw()		
		matLaw:set_elasticity_tensor_orthotropic2d(1e-7, 0, 1, 1)
		elemDisc:set_material_law(matLaw)
		elemDisc:set_quad_order(2) --- !!!!!!!!
		print("--------------------------------------")
		print(elemDisc:config_string())
		
		domainDisc:add(elemDisc)
		
		if util.HasParamOption("-isolator") then
		elemDisc = SmallStrainMechanics(varNames, "Isolator")	
		matLaw = HookeLaw()		
		matLaw:set_elasticity_tensor_orthotropic2d(1, 0, 1, 1)
		elemDisc:set_material_law(matLaw)
		elemDisc:set_quad_order(2) --- !!!!!!!!
		print("--------------------------------------")
		print(elemDisc:config_string())
		
		domainDisc:add(elemDisc)
		end
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

	domainDisc:add(dirichletBND)
	-- domainDisc:add(neumannDisc)
	
	-- ? not necessary
	--OrderCuthillMcKee(approxSpace, false);
	return approxSpace, domainDisc
end

function zero_2(x, y, t)
	return 0
end

function zero_3(x, y, z, t)
	return 0
end

function dirichlet_zero_2(x, y, t)
	return true, 0
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

function IsolatorDiffusion2d(x, y, t)
	--if x > -0.51 and x < 0.51 and y < 0.51 and y > -0.51 then
	if x > -0.5 and x < 0.5 and y < 0.5 and y > -0.5 then 
		return amgEps, 0, 0, amgEps
	else
		return 1, 0, 0, 1
	end
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
	elseif type(tab) ~= "userdata" then		
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
	{ "delta", amg:get_delta() },
	{ "tv", util.GetParam("-amgTV", "const") }
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
	
	local amgTestVcycle = true
	--local amgTestWcycle = true
	--local amgTest2cycle = true
	local ret = {}	
	--print("CG, "..name)  			ret.cg = solve_lin(BiCGStab(), amg, "CG("..name..")")
	if amgTestVcycle then
		amg:set_cycle_type(1)
		print("LinearSolver_V_, "..name)	ret.linV = solve_lin(LinearSolver(), amg, "LS_V("..name..")")
	end
	if amgTestWcycle then
		amg:set_cycle_type(2)
		print("LinearSolver_W_, "..name)	ret.linW = solve_lin(LinearSolver(), amg, "LS_W("..name..")")
	end
	if amgTest2cycle then
		amg:set_Y_cycle(100, 1e-12, 1e-12)	
		print("LinearSolver_2_, "..name)	ret.linY = solve_lin(LinearSolver(), amg, "LS_2("..name..")")
	end
	return ret
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


function test_smoothsteps(smoother, steps)
	local res = {}
	for i, v in pairs(steps) do
		res[v] = solve_amg(smoother, v, smoother:name().."x"..tostring(v))
	end
	return res
end

function AssembleAxb(approxSpace, domainDisc)
	local A = MatrixOperator()
	local b = GridFunction(approxSpace)
	local u = GridFunction(approxSpace)
	u:set(0.0)
	domainDisc:adjust_solution(u)
	
	-- 3. assemble matrix and rhs
	domainDisc:assemble_linear(A, b)
	if util.HasParamOption("-writeMatrices") then
		SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
	end
	return A, u, b
end


function amgsetup.CreateApproxSpaceAndDisc()
	if util.HasParamOption("-LE") then
		local le=nil
		if util.HasParamOption("-isolator") then
			le={Inner={nu=0.4,E=20}, Isolator={nu=0.0, E=20}}
		end
		approxSpace, domainDisc = CreateLEDisc({Domain = dom, BoundaryName = boundaries[1], le=le})
	else
		approxSpace = CreateStdApproxSpace(dom)	
		local diffusion = nil
		local dim = approxSpace:get_dim()
		amgEps = util.GetParamNumber("-amgEps", 1e-7)
		if util.HasParamOption("-isolator") then
			diffusion = {Inner=1, Isolator=amgEps}
		elseif util.HasParamOption("-isolator2") then
			diffusion = {Inner="IsolatorDiffusion"..dim.."d"}
		end
		domainDisc = CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries, diffusion=diffusion })
		
	end
end

function amgsetup.LoadStdDomain2d()
	requiredSubsets = {"Inner", "Boundary"}
	
	if util.HasParamOption("-cembalo") then
		gridName = util.GetParam("-grid", "cembalo/cembalo_platte_2d_sci_1.ugx")
		requiredSubsets = {"Inner", "Boundary"}
		boundaries = {"Boundary"}
	elseif util.HasParamOption("-33k") then
		gridName = util.GetParam("-grid", "amgtests/unstructured_33k.ugx") -- -???
		requiredSubsets = {"Inner", "Left", "Right", "Top", "Bottom"}
		boundaries = {"Bottom"}
	elseif util.HasParamOption("-3k") then
		gridName = util.GetParam("-grid", "amgtests/unstructured_3k.ugx")
		requiredSubsets = {"Inner", "Left", "Right", "Top", "Bottom"}
		boundaries = {"Bottom"}
	elseif util.HasParamOption("-beam") then
		gridName = util.GetParam("-grid", "beamlong.ugx")
		requiredSubsets = {"Inner", "Boundary"}
		boundaries = {"Boundary"}
	elseif util.HasParamOption("-isolator") then
		gridName = "amgtests/unit_square_isolator.ugx"
		requiredSubsets = {"Inner", "Boundary"}
		boundaries = {"Boundary"}
	elseif util.HasParamOption("-squareUG") then
		gridName = "amgtests/unit_square_unstructured.ugx"
		requiredSubsets = {"Inner", "Left", "Right", "Top", "Bottom"}
		boundaries = {"Bottom"}
	else
		gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2.ugx")
		requiredSubsets = {"Inner", "Boundary"}
		boundaries = {"Boundary"}
	end
		
	dom = amgsetup.MyLoadGrid({ dim=2, cpu = cpu, gridName = gridName, requiredSubsets = requiredSubsets} )
end

function amgsetup.LoadStdDomain3d()
	requiredSubsets = {"Inner", "Boundary"}
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_2x2x2.ugx") --unit_cube_01_hex_4x4x4.ugx")
	dom = amgsetup.MyLoadGrid({ dim=3, cpu =cpu, gridName = gridName, requiredSubsets = requiredSubsets} )
	boundaries = {"Boundary"}
end



function amgsetup.LoadStdDomAndApproxSpace()
	-- get dim
	dim = 2
	if util.HasParamOption("-3d") then
		dim = 3
	end
	
	-- setup "cpu" field (for blocks)
	if util.HasParamOption("-LE") then		
		cpu = dim
	else
		cpu = 1
	end
	
	if dim == 2 then
		amgsetup.LoadStdDomain2d()
	else
		amgsetup.LoadStdDomain3d()
	end
	
	amgsetup.CreateApproxSpaceAndDisc()
end