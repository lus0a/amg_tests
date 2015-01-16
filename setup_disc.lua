------------------------------------------------------------------------------------------------------------------------------------
-- grid setup
amgdisc = amgdisc or {}
ug_load_script("amgteaching_problems.lua") 

function amgdisc.MyLoadGrid(p)
	InitUG(p.dim, AlgebraType("CPU", p.cpu))	
	return util.EasyLoadGrid(p)
end
------------------------------------------------------------------------------------------------------------------------------------
-- Discretization

function amgdisc.CreateStdDisc(p)
	NeededParameters(p, {"boundaries", "approxSpace"})

	local domainDisc = DomainDiscretization(p.approxSpace)
	local dim = p.approxSpace:get_dim()
	 
	 if p.diffusion == nil then
	 	local elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
		elemDisc:set_diffusion(1.0)
		--elemDisc:set_diffusion("mydiffusion")
		elemDisc:set_source(0) --"util_std_zero_"..dim)
		domainDisc:add(elemDisc)
	else
		for i, v in pairs(p.diffusion) do
			local elemDisc = ConvectionDiffusion("c", i, "fv1")
			elemDisc:set_diffusion(v)
			elemDisc:set_source(0) --"util_std_zero_"..dim)
			domainDisc:add(elemDisc)
		end
	end
		
	local dirichletBND = DirichletBoundary()
	
	bnds = util.GetParamNumber("-bnd", 1)
	for i,v in pairs(p.boundaries) do
		if i <= bnds then
			dirichletBND:add("util_std_dirichlet_zero_"..dim, "c", v)
		end
	end
				
	domainDisc:add(dirichletBND)
	
	return domainDisc
end


function amgdisc.CreateLognormalDisc(p)
	NeededParameters(p, {"boundaries", "approxSpace"})

	local domainDisc = DomainDiscretization(p.approxSpace)
	local dim = p.approxSpace:get_dim()
	
	local rm_nodes = util.GetParamNumber("-RF_N", 100)
	local rm_mean_f = util.GetParamNumber("-RF_mean_f", 0)
	local rm_sigma_f = util.GetParamNumber("-RF_sigma_f", 1)
	local rm_sigma = util.GetParamNumber("-RF_sigma", 0.01)
	local rm_sigma_1 = util.GetParamNumber("-RF_sigma1", rm_sigma)
	local rm_sigma_2 = util.GetParamNumber("-RF_sigma2", rm_sigma)
	
	 
 	local elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
 	srand(GetClockS()*100)
	local NormalField=LognormalRandomField(rm_nodes, rm_mean_f, rm_sigma_f, rm_sigma)
	elemDisc:set_diffusion(NormalField)
	--elemDisc:set_diffusion("mydiffusion")
	elemDisc:set_source(0) --"util_std_zero_"..dim)
	domainDisc:add(elemDisc)
		
	local dirichletBND = DirichletBoundary()
	
	bnds = util.GetParamNumber("-bnd", 1)
	for i,v in pairs(p.boundaries) do
		if i <= bnds then
			dirichletBND:add("util_std_dirichlet_zero_"..dim, "c", v)
		end
	end
				
	domainDisc:add(dirichletBND)
	
	return domainDisc
end

function amgdisc.CreateLEDisc(parameters)
	local dom = parameters.Domain
	local BoundaryName = parameters.BoundaryName
	local dim = dom:get_dim()
	-- choose algebra
	bLinearElasticity = true
	if util.HasParamOption("-laplace") then
		bLinearElasticity = false
	end
	
	InitUG(dim, AlgebraType("CPU", cpu))
	
	
	-- get approxSpace
	if dim == 3 then
		approxSpace = util.EasyCreateApproxSpace(dom, {"ux", "uy", "uz"})
	else
		approxSpace = util.EasyCreateApproxSpace(dom, {"ux", "uy"})
	end
	
	
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
			elemDisc:set_quad_order(3) --- !!!!!!!!
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
			dirichletBND:add("util_std_dirichlet_zero_"..dim, "ux", boundaryName)
			dirichletBND:add("util_std_dirichlet_zero_"..dim, "uy", boundaryName)
			if dim == 3 then
				dirichletBND:add("util_std_dirichlet_zero_"..dim, "uz", boundaryName)
			end
		else
			dirichletBND:add("util_std_dirichlet_zero_"..dim, "c", boundaryName)
		end
	end
	
	bnds = util.GetParamNumber("-bnd", 1)
	for i,v in pairs(parameters.boundaries) do
		if i <= bnds then
			setBoundary(v)
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




function amgdisc.AssembleAxb(approxSpace, domainDisc)
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

function LoadAMGTeaching(approxSpace, boundaries)
	myproblem = util.GetParamFromList("-amgTeaching", "?", amgTeaching.problem)
	print(myproblem)
	elemDisc = amgTeaching.CreateProblemDependentElemDisc(myproblem, 2)
	return amgTeaching.SetUpProblem(approxSpace, myproblem, elemDisc, boundaries)
end


function amgdisc.CreateApproxSpaceAndDisc()
	if util.HasParamOption("-LE") then
		local le=nil
		if util.HasParamOption("-isolator") then
			le={Inner={nu=0.4,E=20}, Isolator={nu=0.0, E=20}}
		end
		approxSpace, domainDisc = amgdisc.CreateLEDisc({Domain = dom, boundaries=boundaries, le=le})
	else
		approxSpace = util.EasyCreateApproxSpace(dom, {"c"})	
		local diffusion = nil
		local dim = approxSpace:get_dim()
		amgEps = util.GetParamNumber("-amgEps", 1e-7)
		if util.HasParamOption("-isolator") then
			diffusion = {Inner=1, Isolator=amgEps}
		elseif util.HasParamOption("-isolator2") then
			diffusion = {Inner="IsolatorDiffusion"..dim.."d"}
		elseif util.GetParam("-amgTeaching", "no") ~= "no" then
			domainDisc = LoadAMGTeaching(approxSpace, boundaries)
			return
		end
		if util.HasParamOption("-lognormal") then
			domainDisc = amgdisc.CreateLognormalDisc({approxSpace=approxSpace, boundaries=boundaries })
		else		
			domainDisc = amgdisc.CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries, diffusion=diffusion })
		end
		
	end
end

function amgdisc.LoadStdDomain2d()
	requiredSubsets = {"Inner", "Boundary"}
	
	local requiredSubsets4 = {"Inner", "Bottom", "Top", "Left", "Right"}
	
	-- select number of used boundaries with -bnd
	-- e.g. -bnd 1 -> only "Bottom", -bnd 4 -> all.
	local boundaries4 = {"Bottom", "Top", "Left", "Right"}
	
	local requiredSubsets2 = {"Inner", "Boundary"}
	local boundaries2 = {"Boundary"}
	
	if util.HasParamOption("-cembalo") then
		gridName = util.GetParam("-grid", "cembalo/cembalo_platte_2d_sci_1.ugx")
		requiredSubsets = requiredSubsets2
		boundaries = boundaries2
	elseif util.HasParamOption("-33k") then
		gridName = util.GetParam("-grid", "amgtests/unstructured_33k.ugx") -- -???
		requiredSubsets = requiredSubsets4
		boundaries = boundaries4		
	elseif util.HasParamOption("-gitarre") then
		gridName = util.GetParam("-grid", "gitarre/gitarre2d.ugx") -- -???
		requiredSubsets = {"Inner", "OuterZargen", "OuterSchallloch"}
		boundaries = {"OuterZargen", "OuterSchallloch"}
	elseif util.HasParamOption("-3k") then
		gridName = util.GetParam("-grid", "amgtests/unstructured_3k.ugx")
		requiredSubsets = requiredSubsets4
		boundaries = boundaries4
	elseif util.HasParamOption("-beam") then
		gridName = util.GetParam("-grid", "beamlong.ugx")
		requiredSubsets = requiredSubsets2
		boundaries = boundaries2
	elseif util.HasParamOption("-isolator") then
		gridName = "amgtests/unit_square_isolator.ugx"
		requiredSubsets = requiredSubsets2
		boundaries = boundaries2
	elseif util.HasParamOption("-squareUG") then
		gridName = "amgtests/unit_square_unstructured.ugx"
		requiredSubsets = requiredSubsets4
		boundaries = boundaries4
	else
		gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2.ugx")
		requiredSubsets = requiredSubsets2
		boundaries = boundaries2
	end
		
	dom = amgdisc.MyLoadGrid({ dim=2, cpu = cpu, gridName = gridName, requiredSubsets = requiredSubsets} )
end

function amgdisc.LoadStdDomain3d()
	requiredSubsets = {"Inner", "Boundary"}
	gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_2x2x2.ugx") --unit_cube_01_hex_4x4x4.ugx")
	dom = amgdisc.MyLoadGrid({ dim=3, cpu =cpu, gridName = gridName, requiredSubsets = requiredSubsets} )
	boundaries = {"Boundary"}
end



function amgdisc.LoadStdDomAndApproxSpace()
	-- get dim
	dim = 2
	if util.HasParamOption("-3d") then
		dim = 3
	end
	
	-- setup "cpu" field (for blocks)
	if util.HasParamOption("-LE") then		
		if util.HasParamOption("-LE1") then		
			cpu = 1
		else
			cpu = dim
		end	
	else
		cpu = 1
	end
	
	if dim == 2 then
		amgdisc.LoadStdDomain2d()
	else
		amgdisc.LoadStdDomain3d()
	end
	
	amgdisc.CreateApproxSpaceAndDisc()
end
----------------------------------------
-- callback functions for disc



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


------------------------------------------------------------------------------------------------------------------------------------
