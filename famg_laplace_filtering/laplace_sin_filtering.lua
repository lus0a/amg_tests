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

gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")

util.CheckAndPrintHelp("");

matrixWritePath=CreateOutputDirectory("RAMmatrices")

dom = MyLoadGrid({ gridName = gridName} )

approxSpace = CreateStdApproxSpace(dom)
domainDisc = CreateStdDisc(approxSpace)



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
domainDisc:assemble_linear(A, b)
SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
ev = {}
function myEVFunction(x, y, t)
	if nuX==0 or nuY==0 then return (1+math.sin(math.pi*nuX*x))*(1+math.sin(math.pi*y*nuY)) end
	return math.sin(math.pi*nuX*x)*math.sin(math.pi*y*nuY)
end




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

amg = FAMGPreconditioner()
SetFAMGCommandLineParameters(amg, {
	--delta=0.030473,
	delta=0.1,
	smooths=4,
	maxLevels=2, maxNodesForBase=1	
	})
amg:set_deficit_solve_tolerance(1e-15)

jac = Jacobi(0.6)
amg:set_presmoother(jac)
amg:set_postsmoother(jac)

amg:set_matrix_write_path(matrixWritePath)
amg:set_write_testvectors(true)
amg:set_delta(0.0303473)
--amg:set_delta(0.03035)
amg:set_delta(0.1)

amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
amg:set_testvector_smooths(4)
amg:set_testvector_smoother(Jacobi(0.6))


amg:set_write_f_values(true)
--amg:set_debug_write_all_vectors(true)
amg:set_fsmoothing(false)
--amg:set_injection_restriction(true)

print(amg:config_string())

-- create Linear Solver
linSolver = LinearSolver(amg, convCheck)
-- choose some solver
print(amg:config_string())

solver = linSolver
	
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

function solve_random()
	local f = GridFunction(approxSpace)
	f:set_random(-1, 1)		
	solver:apply_return_defect(f, ZeroFunction())
	local avgRate = convCheck:avg_rate()
	if avgRate < 1e-12 then
		convMsg = avgRate.." - EXACT!"
	else
		convMsg = avgRate
	end	
	lastReduction = convCheck:defect()/convCheck:previous_defect()
	print("solve random avgConv = "..convMsg.." lastReduction = ".. lastReduction )
	if false and avgRate < 1e-12 then return "EXACT", "EXACT"
	else 
	return string.format("%2.2e", avgRate), string.format("%2.2e", lastReduction) 
	end
end

checkEVs = util.GetParamNumber("-checkEVs", 5)