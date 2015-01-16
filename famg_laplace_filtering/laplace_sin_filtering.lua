-------------------------------------------------------------------------------
--[[!
-- \file amg/laplace_sin_filtering.lua
-- \{
-- We set a Laplace-problem up here
--
]]--
--------------------------------------------------------------------------------


--PrintBuildConfiguration()

ug_load_script("../setup.lua")

InitSignals()
amgdisc.LoadStdDomAndApproxSpace()
A, u, b = amgdisc.AssembleAxb(approxSpace, domainDisc)

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------


convCheck = ConvCheck(30, 1e-11, 1e-12, false)
precond = Jacobi()


--------------------------------------------------------------------------------
-- Solve
--------------------------------------------------------------------------------

-- 1. Create Algebra


u:set(0.0)
domainDisc:adjust_solution(u)

-- 3. assemble matrix and rhs
domainDisc:assemble_linear(A, b)
if util.HasParamOption("-writeMatrices") then
SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
end
ev = {}




evXs = util.GetParamNumber("-evXs", 1)
evYs = util.GetParamNumber("-evYs", 1)
evXstart = util.GetParamNumber("-evXstart", 1)
evYstart = util.GetParamNumber("-evYstart", 1)

tv =  TestvectorsExternal()
for x=evXstart,evXs do
	ev[x] = {}
	for y=evYstart,evYs do
		ev[x][y] = CreateSinGridFunction(x, y)
		tv:add(ev[x][y])
	end
end

function CreateTV(tab)
	local tve = TestvectorsExternal()	
	for i,v in pairs(tab) do
		print("Adding "..v[1]..", "..v[2])
		tve:add(CreateSinGridFunction(v[1], v[2]))
	end
	return tve
end 

function CreateSinVectors(tab)
	local vt = {}	
	for i,v in pairs(tab) do
		print("Adding "..v[1]..", "..v[2])
		table.insert(vt, CreateSinGridFunction(v[1], v[2]) )
	end
	return vt
end 



--amg:set_testvectors(tv)
--amg:set_testvectors( CreateTV({ {1, 1}, {1, 2}, {2,1} }) )
function createAMG1()
	local amg = FAMGPreconditioner()
	amgsetup.SetFAMGCommandLineParameters(amg, {
		--delta=0.030473,
		delta=0.03,
		smooths=4,
		maxLevels=2, maxNodesForBase=1	
		})
	amg:set_deficit_solve_tolerance(1e-15)
	
	jac = Jacobi(0.6)
	amg:set_presmoother(jac)
	amg:set_postsmoother(jac)
	
	if util.HasParamOption("-writeMatrices") then
		matrixWritePath=CreateOutputDirectory("RAMmatrices")
		amg:set_matrix_write_path(matrixWritePath)
	end
	amg:set_write_testvectors(true)
	
	
	amg:set_write_f_values(true)
	--amg:set_debug_write_all_vectors(true)
	amg:set_fsmoothing(false)
	--amg:set_injection_restriction(true)
	
	return amg
end
	
function createLS(amg)
	-- create Linear Solver
	linSolver = LinearSolver(amg, convCheck)
	-- choose some solver
	print(amg:config_string())
	
	solver = linSolver
end
--------------------------------------------------------------------------------------------

function solve_sin(x, y)
	local f = CreateSinGridFunction(x, y)
	domainDisc:adjust_solution(f)
	solver:apply_return_defect(f, b:clone())
	local avgRate = convCheck:avg_rate()
	if avgRate < 1e-12 then
		convMsg = avgRate.." - EXACT!"
	else
		convMsg = avgRate
	end	
	lastReduction = convCheck:defect()/convCheck:previous_defect()
	print("solve nuX = "..x..", nuY = "..nuY.." avgConv = "..convMsg.." lastReduction = ".. lastReduction )
	if false and avgRate < 1e-12 then return "EXACT", "EXACT"
	else 
	return string.format("%2.2e", avgRate), string.format("%2.2e", lastReduction) 
	end
end

checkEVs = util.GetParamNumber("-checkEVs", 5)

--local a = CreateSinGridFunction(1, 0)
--SaveVectorForConnectionViewer(a, "10.vec")
