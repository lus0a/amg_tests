--   Author: 	Arne Naegel
amgTeaching = amgTeaching or {}

-- parameters for diffusion coefficients
epsx = util.GetParamNumber("-epsx", 1.0)
epsy = util.GetParamNumber("-epsy", 1.0)

uvel = util.GetParamNumber("-uvel", 0.0)
vvel = util.GetParamNumber("-vvel", 0.0)


alpha = util.GetParamNumber("-alpha", 0.0)


--------------------------------
-- User Data Functions 2d (begin)
--------------------------------
	
	
-- A) DIFFUSION
-- standard (anisotropic) diffusion 
function stdDiffTensor2d(x, y, t)
--print (epsx.." "..epsy)
return	epsx, 0, 0, epsy end

-- anisotropic diffusion in corners
function anisoDiffTensor2d(x, y, t)
	local fac1 = 1/1+exp(alpha*(x+y-1)) 
	local fac2 = 1- fac1
	return	fac1*epsx+fac2, 0, 
			0, fac1+fac2*epsy
end

-- hedgehog diffusion 
function hedgehogDiffTensor2d(x, y, t)
	if x<0.5 then 
		if y<0.5 then return 1.0, 0, 0, 1.0 
		else return epsx, 0, 0, 1.0 end
	else
		if y<0.5 then return 1.0, 0, 0, epsy 
		else return 1.0, -1.0, -1.0, 1.0 end
	end
	
	-- should never happen
	return 0, 0, 0, 0 
end

	
	
-- discontinuous coefficient
function jumpDiffTensor2d(x, y, t)
	if (math.abs (x-0.5)<1/4 and math.abs(y-0.5)<1/4) then
		return	epsx, 0, 0, epsy
	else 
		return	1, 0, 0, 1
	end
end

-- B) VELOCITIES
function zeroVelocityField2d(x, y, t)
	return	0, 0
end

function constVelocityField2d(x, y, t)
	return	uvel, vvel
end

function recirculatingVelocityField2d(x, y, t)
	local urecirc = 4.0*x*(1-x)*(1-2*y) 
	local vrecirc = -4.0*y*(1-y)*(1-2*x)
	return	urecirc, vrecirc
end

-- C) REACTIONS
function zeroReaction2d(x, y, t)
	return	0
end

	
	-- D) SOURCES 
function stdSource2d(x, y, t)
	local s = 2*math.pi
	return	s*s*(math.sin(s*x) + math.sin(s*y))
	--return -2*y
--		return 0;
end

-- E) BOUNDARY
function stdNeumannBnd2d(x, y, t)
	--local s = 2*math.pi
	--return -s*math.cos(s*x)
	return true, -2*x*y
end

function stdDirichletBnd2d(x, y, t)
	--	local s = 2*math.pi
	--	return true, math.sin(s*x) + math.sin(s*y)
	return true, 0
	--return true, x*x*y
--		return true, 2.5
end
--------------------------------
-- User Data Functions 2d (end)
--------------------------------


--------------------------------------------------------------------------------
--  Setup problems
--------------------------------------------------------------------------------

amgTeaching.problem = {}

-- depending on the dimension we're choosing the appropriate callbacks.
-- we're using the .. operator to assemble the names (dim = 2 -> "ourDiffTensor2d")

-- Problem 1: Poisson's eqn. (w/ anisotropy)
amgTeaching.problem["Poisson"] =
{ 
	desc = "Poisson's eqn. (w/ anisotropy)",
	diffusion = "stdDiffTensor2d",
	velocity = "zeroVelocityField2d",
	reaction = "zeroReaction2d",
	source = "stdSource2d",
	neumann = "stdNeumannBnd2d",
	dirichlet = "stdDirichletBnd2d",
	configs=
	{ 
		{["epsx"]= 1.0, ["epsy"]= 1.0},
		{["epsx"]= 1.0, ["epsy"]= 1e-1},
		{["epsx"]= 1.0, ["epsy"]= 0.075},
		{["epsx"]= 1.0, ["epsy"]= 5e-2},
		{["epsx"]= 1.0, ["epsy"]= 2e-2},
		{["epsx"]= 1.0, ["epsy"]= 1e-2},
		{["epsx"]= 1.0, ["epsy"]= 5e-3},
		{["epsx"]= 1.0, ["epsy"]= 1e-3},
		{["epsx"]= 1.0, ["epsy"]= 5e-4},
		{["epsx"]= 1.0, ["epsy"]= 1e-4},
		{["epsx"]= 1.0, ["epsy"]= 1e-5},
		{["epsx"]= 1.0, ["epsy"]= 1e-6},
	--{["epsy"]= 1.0, ["epsx"]= 1e-1},
	--{["epsy"]= 1.0, ["epsx"]= 1e-2},
	--{["epsy"]= 1.0, ["epsx"]= 1e-3},
	--{["epsy"]= 1.0, ["epsx"]= 1e-4},
	--{["epsy"]= 1.0, ["epsx"]= 1e-5},
	--{["epsy"]= 1.0, ["epsx"]= 1e-6}
	}
}

-- Problem 2: Convection-diffusion (w/ anisotropy) 
amgTeaching.problem["ConvDiff"] = 
{
	desc = "Convection-diffusion (w/ anisotropy)",
	diffusion = "stdDiffTensor2d",
	velocity = "constVelocityField2d",
	reaction = "zeroReaction2d",
	source = "stdSource2d",
	neumann = "stdNeumannBnd2d",
	dirichlet = "stdDirichletBnd2d",
	configs={
		{["epsx"]= 1.0, ["epsy"]= 1.0},
		{["epsx"]= 1e-1, ["epsy"]= 1e-1},
		{["epsx"]= 1e-2, ["epsy"]= 1e-2},
		{["epsx"]= 1e-3, ["epsy"]= 1e-3},
		{["epsx"]= 1e-4, ["epsy"]= 1e-4},
		{["epsx"]= 1e-5, ["epsy"]= 1e-5},
		{["epsx"]= 1e-6, ["epsy"]= 1e-6}
	}
}

-- Problem 3: Recirculating flow (w/ anisotropy) 
amgTeaching.problem["Recirc"] =
{
	desc = "Recirculating flow (w/ anisotropy)",
	diffusion = "stdDiffTensor2d",
	velocity = "recirculatingVelocityField2d",
	reaction = "zeroReaction2d",
	source = "stdSource2d",
	neumann = "stdNeumannBnd2d",
	dirichlet = "stdDirichletBnd2d",
	configs ={
		{["epsx"]= 1.0, ["epsy"]= 1.0},
		{["epsx"]= 1.0, ["epsy"]= 1.0},
		{["epsx"]= 1e-1, ["epsy"]= 1e-1},
		{["epsx"]= 1e-2, ["epsy"]= 1e-2},
		{["epsx"]= 1e-3, ["epsy"]= 1e-3},
		{["epsx"]= 1e-4, ["epsy"]= 1e-4},
		{["epsx"]= 1e-5, ["epsy"]= 1e-5},
		{["epsx"]= 1e-6, ["epsy"]= 1e-6}
	}
}

-- Problem 4: Hedgehog diffusion (w/ anisotropy, four different subdomains) 
amgTeaching.problem["Hedgehog"] = 
{
	desc = "Hedgehog-type Poisson's eqn. (w/ anisotropy)",
	diffusion = "hedgehogDiffTensor2d",
	velocity = "zeroVelocityField2d",
	reaction = "zeroReaction2d",
	source = "stdSource2d",
	neumann = "stdNeumannBnd2d",
	dirichlet = "stdDirichletBnd2d",
	configs = {
		{["epsx"]= 1.0, ["epsy"]= 1.0},
		{["epsx"]= 1e-1, ["epsy"]= 1e-1},
		{["epsx"]= 1e-2, ["epsy"]= 1e-2},
		{["epsx"]= 1e-3, ["epsy"]= 1e-3},
		{["epsx"]= 1e-4, ["epsy"]= 1e-4},
		{["epsx"]= 1e-5, ["epsy"]= 1e-5},
		{["epsx"]= 1e-6, ["epsy"]= 1e-6}
	}
}

-- Problem 5: Isolator problem (embeded squarent subdomains) 
amgTeaching.problem["Isolator"] =
{
	desc = "Isolator problem (jumping coefficient)",
	diffusion = "jumpDiffTensor2d",
	velocity = "zeroVelocityField2d",
	reaction = "zeroReaction2d",
	source = "stdSource2d",
	neumann = "stdNeumannBnd2d",
	dirichlet ="stdDirichletBnd2d",
	configs = 
	{
		{["epsx"]= 1e+6, ["epsy"]= 1e+6},
		{["epsx"]= 1e+4, ["epsy"]= 1e+4},
		{["epsx"]= 1e+2, ["epsy"]= 1e+2},
		{["epsx"]= 1e+1, ["epsy"]= 1e+1},
		{["epsx"]= 1.0, ["epsy"]= 1.0},
		{["epsx"]= 1e-1, ["epsy"]= 1e-1},
		{["epsx"]= 1e-2, ["epsy"]= 1e-2},
		{["epsx"]= 1e-4, ["epsy"]= 1e-4},
		{["epsx"]= 1e-6, ["epsy"]= 1e-6}
	}
}

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------	

function amgTeaching.CreateProblemDependentElemDisc(problem, dim)
	print ("Setting up Assembling")
	
	-- Select upwind
	if dim == 2 then 
	
	--upwind = NoUpwind2d()
	upwind = FullUpwind2d()
	--upwind = WeightedUpwind2d(); upwind:set_weight(0.0)
	--upwind = PartialUpwind2d()
	else print("Dim not supported for upwind"); exit() end
	
	elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
	elemDisc:set_upwind(upwind)
	elemDisc:set_diffusion(LuaUserMatrix(problem.diffusion))
	elemDisc:set_velocity(LuaUserVector(problem.velocity))
	elemDisc:set_reaction_rate(LuaUserNumber(problem.reaction))
	elemDisc:set_source(LuaUserNumber(problem.source))
	return elemDisc
end

function amgTeaching.SetUpProblem(approxSpace, myproblem, elemDisc, boundaries)
	-----------------------------------------------------------------
	--  Setup Neumann Boundary
	-----------------------------------------------------------------
	
	--neumannDisc = NeumannBoundary("c")
	--neumannDisc:add(neumann, "NeumannBoundary", "Inner")
	
	-----------------------------------------------------------------
	--  Setup Dirichlet Boundary
	-----------------------------------------------------------------
	
	if boundaries == nil then
		boundaries = {"DirichletBoundary"}
	end
	
	dirichletBND = DirichletBoundary()
	for k, v in pairs(boundaries) do
		dirichletBND:add(LuaCondUserNumber(myproblem.dirichlet), "c", v)
	end
	
	-------------------------------------------
	--  Setup Domain Discretization
	-------------------------------------------
	
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(elemDisc)
	--domainDisc:add(neumannDisc)
	domainDisc:add(dirichletBND)
	
	
	return domainDisc
end

function amgTeaching.AssembleAxb(approxSpace, domainDisc)
-- create operator from discretization
	local A = AssembledLinearOperator()
	A:set_discretization(domainDisc)
	
	-- get grid function
	u = GridFunction(approxSpace)
	b = GridFunction(approxSpace)
	
	-- 0. Reset start solution
	u:set(0.0)
	-- 1. init operator
	--print("Init operator (i.e. assemble matrix).")
	if AssembleLinearOperatorRhsAndSolution(A, u, b) == false then 
		print("Could not assemble operator"); exit(); 
	end
	return A, x, b
end

