-------------------------------------------------------------------------------
--[[!
-- \file amg/cycle_check.lua
-- this script tests the cycle types.
-- main function is check_cycle(k, eta, gamma)
-- k: is the number of multigrid levels
-- eta: is the convergence rate of the base solver on level k
--     (i.e. Jacobi iteration solving up to reduction of eta)
-- gamma: cycle type.
--
]]--
--------------------------------------------------------------------------------



SetMinSecondsUntilProgress(3)
--PrintBuildConfiguration()

ug_load_script("setup.lua")

GetLogAssistant():enable_file_output(true, "runinfo.txt")
print("COMMAND LINE: ")
print(util.GetCommandLine())

InitSignals()

outfile = util.GetParam("-outfile", nil, "output directory")

maxLevels = util.GetParamNumber("-maxLevels", 2)
numSmooths = util.GetParamNumber("-numSmooths", 2)

if true then
	amgdisc.LoadStdDomAndApproxSpace()
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



--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

precond = Jacobi()
OrderCuthillMcKee(approxSpace, true)

--------------------------------------------------------------------------------
-- Solve
--------------------------------------------------------------------------------

-- 1. Create Algebra
A, u, b = amgdisc.AssembleAxb(approxSpace, domainDisc)

function getGMG()
	local amg = RSAMGPreconditioner()
	return amg
end

function getMyAMG(k, gamma)

if false then
	amg = RSAMGPreconditioner()
else
	amg = FAMGPreconditioner()
	amgsetup.SetFAMGCommandLineParameters(amg, {
		smooths=2,
		--smoother="jac",
		smoother="jac",
		maxNodesForBase=2,
		maxLevels=2,
		delta=0.01	
		}) 
	amg:set_write_f_values(true)
end

amg:set_max_levels(k)
amg:set_cycle_type(gamma)
amg:set_max_nodes_for_base(2)
--amg:set_Y_cycle(4, 1e-12, 0.07)
	
return amg
	
end
	

function get_two_grid_correction_rate(maxLevels)
	--[[ return 0.32
end
function blaaa()
]]--
	amg = getMyAMG(maxLevels, 1)
	
	gmgSmoother = Jacobi(0.5)
	amg:set_postsmoother(gmgSmoother)
	amg:set_presmoother(gmgSmoother)
	--if util.HasParamOption("-writeMatrices") then
		matrixWritePath=CreateOutputDirectory("RAMmatrices")
		
		uuu = GridFunction(approxSpace)
		uuuu = GridFunctionPositionProvider(uuu)
		amg:set_position_provider( uuuu )
		--matrixWritePath = "/Users/mrupp/matrices/"
		amg:set_matrix_write_path(matrixWritePath)
		--amg:set_write_f_values(true)
		--amg:set_debug_write_all_vectors(true)
	--end

	
	gmg = amg
	
	gmg:set_num_presmooth(numSmooths)
	gmg:set_num_postsmooth(numSmooths)
	
	gmg:set_base_solver(LU())
	--local bcc = ConvCheck(50000, 1e-20, 0.1, false)
	--bcc:set_supress_unsuccessful(false)
	--gmg:set_base_solver(LinearSolver(Jacobi(), bcc))
	
	local gmgConvCheck = ConvCheck(100, 1e-12, 1e-12, true)
	local linSolverGMG = LinearSolver(gmg, gmgConvCheck)
	linSolverGMG:init(A)
	gmgRes=solve_random_d(linSolverGMG, approxSpace, gmgConvCheck)
	--gmgRes=solve_random_d(linSolverGMG, approxSpace, gmgConvCheck)
	
	print(amg:config_string())
	
		amg:set_preiterations_mimum_defect_at_check(1e-30)	
	amg:set_nr_of_preiterations_at_check(20)
	amg:check()
	gmgRes=solve_random_d(linSolverGMG, approxSpace, gmgConvCheck)
	return gmgRes.lastReduction
end

gmgTwoGridConvergence = get_two_grid_correction_rate(maxLevels)
print("Estimated Two-Grid-Convergence = "..gmgTwoGridConvergence)

results = StringTable()
resi=1
function check_cycle(k, eta, gamma)
	if gamma == "F" then return end
	amg = getMyAMG(k+1, gamma)
	
	gmgSmoother = Jacobi(0.5)
	amg:set_postsmoother(gmgSmoother)
	amg:set_presmoother(gmgSmoother)
	
	gmg = amg
	
	gmg:set_num_presmooth(numSmooths)
	gmg:set_num_postsmooth(numSmooths)
	
	local bcc = ConvCheck(50000, 1e-20, eta, false)
	bcc:set_supress_unsuccessful(false)
	--gmg:set_base_solver(LinearSolver(Jacobi(), bcc))
	--gmg:set_base_solver(LinearSolver(GaussSeidel(), bcc))
	gmg:set_base_solver(LinearSolver(Jacobi(), bcc))
	
	
	print(" k = "..k..", eta = "..eta..", gamma = "..gamma)
	print("------------------------------------------------------------")
	
	if gamma == "F" then
		estRed = math.pow(eta, k+1)		
	else
		estRed = math.pow(eta, math.pow(gamma,k))	
		print("gamma^k = "..math.pow(gamma,k)..", eta^(gamma^k) = "..estRed)
	end
	
	if estRed < gmgTwoGridConvergence then
		estRed = gmgTwoGridConvergence
		print("est. reduction = two-grid-reduction")
		star = "*"
	else
		print("est. reduction = coarse-grid-reduction")
		star = ""
	end
	
	
	--gmg = Jacobi(0.66)
	local gmgConvCheck = ConvCheck(100, 1e-12, 1e-12, false)
	local linSolverGMG = LinearSolver(gmg, gmgConvCheck)
	linSolverGMG:init(A)
	gmgRes=solve_random_d(linSolverGMG, approxSpace, gmgConvCheck)
	--print(gmgRes)
	
	print("\nlastReduction        = "..gmgRes.lastReduction)
	print("calculated Reduction = "..estRed)
	
	print("")
	results:set(0, 0, "$k$")
	results:set(0, 1, "$\\eta_H$")
	results:set(0, 2, "$\\gamma$")
	results:set(0, 3, "KonvRate")
	results:set(0, 4, "Abschätzung")
	results:set(0, 5, "$\\eta_H^(\\gamma^k)$")
	
	strReduction = string.format("%5.4f", gmgRes.lastReduction)
	strEstReduction = string.format("%5.4f", estRed)
	
	
	--f:write("k#eta#gamma#reduction#calculatedReduction\n")
	f:write(k.."\t"..eta.."\t"..gamma.."\t"..strReduction.."\t"..strEstReduction.."\n")
	
	results:set(resi, 0, k)
	results:set(resi, 1, eta)
	results:set(resi, 2, gamma)
	results:set(resi, 3, strReduction)
	results:set(resi, 4, strEstReduction..star)
	if gamma=="F" then
		results:set(resi, 5, string.format("%5.4f", math.pow(eta, k+1)))
	else
		results:set(resi, 5, string.format("%5.4f", math.pow(eta, math.pow(gamma,k))))
	end
	resi = resi+1
end


-- !!!
-- für k = 0 muss
-- check_cycle(k=0, eta, gamma) = eta sein, 
-- da dann direkt der basesolver verwendet wird
-- das feinste Level ist also Level 0, k = 1 -> 2-Gitter, k=2 -> 3-Gitter usw. 

if false then
	check_cycle(1, 0.4, 1)
	check_cycle(2, 0.4, 1)
	check_cycle(3, 0.4, 1)
	
	
	check_cycle(1, 0.5, 1)
	check_cycle(2, 0.5, 1)
	check_cycle(3, 0.5, 1)
	
	check_cycle(1, 0.5, "F")
	check_cycle(2, 0.5, "F")
	
	
	check_cycle(1, 0.6, 1)
	check_cycle(2, 0.6, 1)
	check_cycle(3, 0.6, 1)
	
	check_cycle(1, 0.6, 2)
	check_cycle(2, 0.6, 2)
	
	check_cycle(1, 0.6, "F")
	check_cycle(2, 0.6, "F")
	check_cycle(3, 0.6, "F")
	
	check_cycle(1, 0.7, 3)
	
	check_cycle(1, 0.8, 2)
	check_cycle(2, 0.8, 2)
	check_cycle(3, 0.8, 2)
	
	check_cycle(1, 0.8, 3)
	
	check_cycle(1, 0.8, "F")
	check_cycle(2, 0.8, "F")
	check_cycle(3, 0.8, "F")
	
	check_cycle(1, 0.9, 2)
	check_cycle(2, 0.9, 2)
	check_cycle(3, 0.9, 2)
	
	check_cycle(1, 0.9, 3)
	check_cycle(2, 0.9, 3)
	
	check_cycle(1, 0.9, "F")
	check_cycle(2, 0.9, "F")
	check_cycle(3, 0.9, "F")
	check_cycle(4, 0.9, "F")
end

function check_for(iStart, iEnd, iStep, k, gamma, numSteps)
	if gamma == "F" then
		estStart = math.exp(math.log(gmgTwoGridConvergence)/(k+1))
	else
		estStart = math.exp(math.log(gmgTwoGridConvergence)/(math.pow(gamma, k)))
	end
	
	
	--startEta = math.max(estStart-iStep*2, iStart)
	if estStart > 0 then		
		iStart = math.max(estStart-4*iStep, 0.01)
		iStep = (iEnd-iStart)/(numSteps+4)
	else
		iStep = (iEnd-iStart)/(numSteps+4)
	end
	print("Estimated Start is at ".. estStart)
	print("iStart = ".. iStart..", iEnd = "..iEnd)
	print("step = ".. iStep.." numSteps = "..numSteps)
	f = io.open("cycle_check_k_"..k.."_gamma_"..gamma..".data", "w")
	f:write("#eta\treduction\tcalculatedReduction\n")
	for i=iStart,iEnd,iStep do check_cycle(k, i, gamma) end
end


if maxLevels == 2 then
-- 7
	check_for(0.1, 0.99, 0.05, 1, 1, 10)
	check_for(0.1, 0.99, 0.05, 1, 2, 10)
	check_for(0.5, 0.99, 0.02, 1, 3, 10)
	--check_for(0.5, 0.99, 0.05, 1, "F", 10)
end

if maxLevels == 3 then
-- 8
	check_for(0.1, 0.99, 0.05, 2, 1, 10)
	check_for(0.4, 0.99, 0.05, 2, 2, 10)
	check_for(0.7, 0.99, 0.02, 2, 3, 10)
	--check_for(0.5, 0.99, 0.05, 2, "F", 10)
end

if maxLevels == 4 then
-- numrefs 9
	check_for(0.1, 0.99, 0.05, 3, 1, 10)
	check_for(0.7, 0.99, 0.02, 3, 2, 10)
	check_for(0.97, 0.99, 0.001, 3, 3, 10)
	--check_for(0.5, 0.99, 0.05, 3, "F", 10)	
end



results:set_row_seperator(0, "-")
results:set_row_seperator(1, "-")
results:set_row_seperator(results:num_rows(), "-")

print(results:to_latex())



-- geht: ugshell -ex amg/cycle_check_amg.lua -numRefs 0 -3k -writeMatrices -amgXC -amgTV lin
-- geht nicht: ugshell -ex amg/cycle_check_amg.lua -numRefs 0 -3k -writeMatrices -amgXC -amgTV const !!! wg const TV
--[[
geht:
ugshell -ex amg/cycle_check_amg.lua -numRefs 0 -3k -writeMatrices -amgXC -amgTV lin -maxLevels 2
ugshell -ex amg/cycle_check_amg.lua -numRefs 0 -3k -writeMatrices -amgXC -amgTV lin -maxLevels 3
ugshell -ex amg/cycle_check_amg.lua -numRefs 0 -3k -writeMatrices -amgXC -amgTV lin -maxLevels 4


ugshell -ex amg/cycle_check_amg.lua -numRefs 6 -writeMatrices -amgTV lin -amgXC -maxLevels 4
geht auch mit sgs / ilu als basesolver!!!
next step: mit amg als basesolver.

mit numSmooths 1 andere Ergebnisse, da konv-rate schon von anfang an schlechter, dadurch nicht das verhalten, dass bei TV const und kleinem eta die schätzung nicht stimmt
dafür hinterher nciht mehr so genau


test für conv auf tieferen leveln
ugshell -ex amg/cycle_check_amg.lua -maxLevels 3 -numRefs 6 -amgTV const -amgXC -amgDelta 0.3 
]]--