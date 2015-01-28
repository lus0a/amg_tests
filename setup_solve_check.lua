
function last_random_conv(solver, approxSpace, convCheck)
	return solve_random_d(solver, approxSpace, convCheck).lastReduction
end

function solve_random(solver, approxSpace, convCheck)
	local avgRate
	local lastReduction
	local convMsg
  	local res = solve_random_d(solver, approxSpace, convCheck)
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
	--local lastReduction = 1 --convCheck:rate()
	local timeSolveS
	if tBeginS then
		timeSolveS = GetClockS()-tBeginS
	else
		timeSolveS =nil
	end	
	return {avgRate=avgRate, lastReduction=lastReduction, timeSolveS=timeSolveS, steps=convCheck:step()}
end

function solve_random_d(solver, approxSpace, convCheck)
	local tBeginS = GetClockS()
	local f = GridFunction(approxSpace)
	srand(0)
	f:set_random(-1, 1)		
	local d = ZeroFunction()
	solver:apply_return_defect(f, d)
	local cc = solver:convergence_check()		
	local res= GetConvCheckTable(convCheck, tBeginS)
	res.u = f
	res.d = d
	return res		
end



function solve_lin(linSolver, amg, convCheck, name)
	linSolver:set_convergence_check(convCheck)	
	linSolver:set_preconditioner(amg)
	linSolver:init(A)
	return solve_random_d(linSolver, approxSpace, convCheck)
	
end

function solve_amg(smoother, numSmooths, convCheck, name)


	amg:set_presmoother(smoother)
	amg:set_postsmoother(smoother)
	amg:set_num_presmooth(numSmooths)
	amg:set_num_postsmooth(numSmooths)
	
	cgSolver = BiCGStab()
	cgSolver:set_preconditioner(amg)	
	amg:init_smoothers()
	
	print(name..": "..numSmooths.."x "..smoother:config_string())
	
	--local amgTestVcycle = true
	--local amgTestWcycle = true
	--local amgTest2cycle = true
	local ret = {}	
	--print("CG, "..name)  			ret.cg = solve_lin(BiCGStab(), amg, "CG("..name..")")
	if amgTestVcycle then
		amg:set_cycle_type(1)		
		print("LinearSolver_V_, "..name)	ret.linV = solve_lin(LinearSolver(), amg, convCheck, "LS_V("..name..")")
	end
	if amgTestWcycle then
		amg:set_cycle_type(2)
		print("LinearSolver_W_, "..name)	ret.linW = solve_lin(LinearSolver(), amg, convCheck, "LS_W("..name..")")
	end
	if amgTest2cycle then
		--amg:set_Y_cycle(100, 1e-20, 1e-5)
		amg:clear_level_conv_rates()
		amg:set_Y_cycle(100, 1e-20, 0.2)			
		print("LinearSolver_2_, "..name)	
		ret.linY020 = solve_lin(LinearSolver(), amg, convCheck, "LS_2("..name..")")
		print(amg:level_conv_rates())
		
		amg:clear_level_conv_rates()
		amg:set_Y_cycle(100, 1e-20, 0.1)			
		print("LinearSolver_2_, "..name)	
		ret.linY010 = solve_lin(LinearSolver(), amg, convCheck, "LS_2("..name..")")
		print(amg:level_conv_rates())
		
		amg:clear_level_conv_rates()
		amg:set_Y_cycle(100, 1e-20, 0.05)			
		print("LinearSolver_2_, "..name)	
		ret.linY005 = solve_lin(LinearSolver(), amg, convCheck, "LS_2("..name..")")
		print(amg:level_conv_rates())
		
		amg:clear_level_conv_rates()
		amg:set_Y_cycle(100, 1e-20, 0.01)			
		print("LinearSolver_2_, "..name)	
		ret.linY001 = solve_lin(LinearSolver(), amg, convCheck, "LS_2("..name..")")
		print(amg:level_conv_rates())
		
		amg:clear_level_conv_rates()
		amg:set_Y_cycle(100, 1e-20, 0.001)			
		print("LinearSolver_2_, "..name)	
		ret.linY0001 = solve_lin(LinearSolver(), amg, convCheck, "LS_2("..name..")")
		print(amg:level_conv_rates())
	end
	amg:set_cycle_type(1)		
	
	return ret
end



function test_smoothsteps(smoother, convCheck, steps)
	local res = {}
	for i, v in pairs(steps) do
		res[v] = solve_amg(smoother, v, convCheck, smoother:name().."x"..tostring(v))
	end
	return res
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