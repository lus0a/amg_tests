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
	local gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(util.GetBestLU())
	
	
	gmgSmoother = Jacobi(0.5)
	gmg:set_smoother(gmgSmoother);
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)
	gmg:set_rap(false)
	return gmg
end

function get_two_grid_correction_rate()
	gmg = getGMG()
	gmg:set_base_level(0)
	gmg:set_cycle_type(3)	
	local gmgConvCheck = ConvCheck(100, 1e-12, 1e-12, false)
	local linSolverGMG = LinearSolver(gmg, gmgConvCheck)
	linSolverGMG:init(A)
	gmgRes=solve_random_d(linSolverGMG, approxSpace, gmgConvCheck)
	return gmgRes.lastReduction
end

gmgTwoGridConvergence = get_two_grid_correction_rate()
print("Estimated Two-Grid-Convergence = "..gmgTwoGridConvergence)

results = StringTable()
resi=1
function check_cycle(k, eta, gamma)
	gmg = getGMG()	
	gmg:set_base_level(numRefs-k-1)
	gmg:set_cycle_type(gamma)
	
	local bcc = ConvCheck(50000, 1e-20, eta, false)
	bcc:set_supress_unsuccessful(false)
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
	
	print("lastReduction        = "..gmgRes.lastReduction)
	print("calculated Reduction = "..estRed)
	
	print("")
	results:set(0, 0, "$k$")
	results:set(0, 1, "$\\eta_H$")
	results:set(0, 2, "$\\gamma$")
	results:set(0, 3, "KonvRate")
	results:set(0, 4, "AbschŠtzung")
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
-- fŸr k = 0 muss
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
	
	print("Estimated Start is at ".. estStart)
	--startEta = math.max(estStart-iStep*2, iStart)
	if estStart > 0 then
		iStep = (iEnd-estStart)/numSteps
		iStart = math.max(estStart-4*iStep, 1e-12)
	else
		iStep = (iEnd-iStart)/(numSteps+4)
	end
	f = io.open("cycle_check_k_"..k.."_gamma_"..gamma..".data", "w")
	f:write("#eta\treduction\tcalculatedReduction\n")
	for i=iStart,iEnd,iStep do check_cycle(k, i, gamma) end
end


if util.HasParamOption("-k1") then
-- 7
	check_for(0.1, 0.99, 0.05, 1, 1, 10)
	check_for(0.1, 0.99, 0.05, 1, 2, 10)
	check_for(0.5, 0.99, 0.02, 1, 3, 10)
	check_for(0.5, 0.99, 0.05, 1, "F", 10)
end

if util.HasParamOption("-k2") then
-- 8
	check_for(0.3, 0.99, 0.05, 2, 1, 10)
	check_for(0.4, 0.99, 0.05, 2, 2, 10)
	check_for(0.7, 0.99, 0.02, 2, 3, 15)
	check_for(0.5, 0.99, 0.05, 2, "F", 10)
end

if util.HasParamOption("-k3") then
-- numrefs 9
	check_for(0.7, 0.99, 0.05, 3, 1, 15)
	check_for(0.7, 0.99, 0.02, 3, 2, 15)
	check_for(0.97, 0.999, 0.001, 3, 3, 20)
	check_for(0.7, 0.99, 0.05, 3, "F", 15)	
end



results:set_row_seperator(0, "-")
results:set_row_seperator(1, "-")
results:set_row_seperator(results:num_rows(), "-")

print(results:to_latex())

--[[
max(x,y) = (x > y ? x : y)
P(x, gamma) = max(x**gamma, 0.191308)
plot "cycle_check_k_1_gamma_1.data" using 2:4, P(x, 1), "cycle_check_k_1_gamma_2.data" using 2:4, P(x, 2), "cycle_check_k_1_gamma_3.data" using 2:4, P(x, 3), "cycle_check_k_1_gamma_F.data" using 2:4, P(x, 2)
plot "cycle_check_k_2_gamma_1.data" using 2:4, P(x, 1), "cycle_check_k_2_gamma_2.data" using 2:4, P(x, 4), "cycle_check_k_2_gamma_3.data" using 2:4, P(x, 9),  "cycle_check_k_2_gamma_F.data" using 2:4, P(x, 3)
plot "cycle_check_k_3_gamma_1.data" using 2:4, P(x, 1), "cycle_check_k_3_gamma_2.data" using 2:4, P(x, 8), "cycle_check_k_3_gamma_3.data" using 2:4, P(x, 27),  "cycle_check_k_3_gamma_F.data" using 2:4, P(x, 4)


plot "cycle_check_k_2_gamma_2.data" using 2:4, P(x, 4)

ugshell -ex amg/cycle_check.lua -numRefs 7 -k1; ugshell -ex amg/cycle_check.lua -numRefs 8 -k2; ugshell -ex amg/cycle_check.lua -numRefs 9 -k3
klein:
ugshell -ex amg/cycle_check.lua -numRefs 6 -k1; ugshell -ex amg/cycle_check.lua -numRefs 7 -k2; ugshell -ex amg/cycle_check.lua -numRefs 8 -k3

]]--

-- auch: 
--[[
 ugshell -ex amg/cycle_check.lua -numRefs 2 -3k -k1
 ugshell -ex amg/cycle_check.lua -numRefs 3 -3k -k2
 ugshell -ex amg/cycle_check.lua -numRefs 4 -3k -k3
 
 
 ugshell -ex amg/cycle_check.lua -numRefs 2 -3k -k1; ugshell -ex amg/cycle_check.lua -numRefs 3 -3k -k2; ugshell -ex amg/cycle_check.lua -numRefs 4 -3k -k3 
 
 ]]--
-- ugshell -ex amg/cycle_check.lua -numRefs 2 -3k



