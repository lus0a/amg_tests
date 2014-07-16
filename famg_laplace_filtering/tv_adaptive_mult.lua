-- hello

ug_load_script("laplace_sin_filtering.lua")
outdir = util.GetParam("-outdir")

checkEVs = util.GetParamNumber("-checkEVs", 8)
function createAMG()
	amg = createAMG1()
	if false then
	amg = FAMGPreconditioner()
	end
	--amg:set_delta(0.031)
	amg:set_delta(0.3)
	amg:set_max_nodes_for_base(2000)
	amg:set_max_levels(20)
	
	numSmooths = 4 --util.GetParamNumber("-amgSmooths", 2)
	
	amg:set_num_presmooth(numSmooths)
	amg:set_num_postsmooth(numSmooths)
	
	smoother = SymmetricGaussSeidel()
	amg:set_presmoother(smoother)
	amg:set_postsmoother(smoother)
	
	---{{{
		--! this section is for {0,0} do be interpolated exactly
		--! NOTE: {0,0} is an exception! it is NOT smooth
		--! (see at dirichlet points)	
		
		--! we MUST set testvectorsmooths = presmooths 
		amg:set_testvector_smooths(amg:get_num_presmooth())
		
		--! of course, we also have to set the same smoother
		amg:set_testvector_smoother(smoother)
		
		--! if we use fsmoothing, the {0,0} testvector will be changed.
		--! however, in the construction of FAMG, we don't know in advance which nodes are fine,
		--! so we can't use f-smoothing here.	
		amg:set_fsmoothing(false)
			
		-- this is no necessary. we can use R = P^T.
		--amg:set_injection_restriction(false)
		
		--! note that all this is not necessary for "real" eigenmodes like sin,
		--! {1,1} for example is so smooth that the smoother hardly does anything	
	---}}}
	
	
	
	--amg:set_testvector_smooths(0)
	amg:set_injection_restriction(false)
	return amg
end


function test(tv, name)
	print("---------------------------------------------------------------")
	print("                        "..name)
	print("---------------------------------------------------------------")
	
	amg:set_testvectors(tv)
	solver:init(A)
	print(amg:level_information_str())
	
	amg:set_num_presmooth(4)
	amg:set_num_postsmooth(4)
	last1 = last_random_conv(solver, approxSpace)
	
	if last1 < 0.3 then
		amg:set_num_presmooth(10)
		amg:set_num_postsmooth(10)
		last = last_random_conv(solver, approxSpace)
		
		if last < 0.2 then
			amg:set_num_presmooth(40)
			amg:set_num_postsmooth(40)
			last = last_random_conv(solver, approxSpace)
			
			if big == false then
				amg:set_num_presmooth(100)
				amg:set_num_postsmooth(100)
				last = last_random_conv(solver, approxSpace)
			end
		end
	end
	return last1
end

function test_ext(tv, name)
print(tv)	
	return test(TestvectorsExternal(tv), name)
end

function test2(tv)
	return test(tv, tv:config_string())
end




if false then

v1 = RandomFunction()
smooth(v1, ILU(), 100, A)
v2 = RandomFunction()
smooth(v2, ILU(), 100, A)
v3 = RandomFunction()
smooth(v3, ILU(), 100, A)
VecScaleAssign(v1, 1/VecNorm(v1), v1)
VecScaleAssign(v2, 1/VecNorm(v2), v2)
VecScaleAssign(v3, 1/VecNorm(v3), v3)
--amg:add_watchnode(0,0)
test_ext({v1}, "pevp")
do return end
end

convCheck:set_verbose(true)

v1 = RandomFunction()
smooth(v1, Jacobi(0.8), 4, A)

lip = LinearIteratorProduct()
lsLIP = LinearSolver(lip, convCheck)

for filteringStep=1,20 do
	amg = createAMG()
	amg:set_one_init(true)
	createLS(amg)
			
	if util.HasParamOption("-writeMatrices") then
		SaveVectorForConnectionViewer(v1, "v"..filteringStep..".vec")
		amg:set_matrix_write_path(CreateOutputDirectory("RAMmatrices/"..filteringStep))
	end
	
	amg:set_testvectors(TestvectorsExternal({v1}))
	solver:init(A)

	lip:add_iterator(amg)
	--aaa=FAMGPreconditioner()
	print(amg:level_information_str())
	
	
	lsLIP:init(A)
	
	if filteringStep ~= 1 then amgRes = solve_random_d(solver, approxSpace) end
	LIPRes = solve_random_d(lsLIP, approxSpace)	
	if LIPRes.lastReduction < 0.1 then break end
	
	v1 = RandomFunction()
	smooth(v1, lip, 10, A)
	VecScaleAssign(v1, 1/VecNorm(v1), v1)	
	
end



if ProcRank() == 0 then
	stats = {
	
	{ "procs", NumProcs() },
	{ "Type", "laplace" },
	{ "numVecs", 1},
	--{ "numPreRefs", numPreRefs},
	{ "numRefs", numRefs },
	{ "dim", dim},
	{ "gridName", gridName},
	{ "main [ms]", GetProfileNode("main"):get_avg_total_time_ms()},
	{ "perform_refinement [ms]", GetProfileNode("perform_refinement"):get_avg_total_time_ms()},
	 
	{ "date", os.date("y%Ym%md%d.%H:%M:%S")},
	{ "SVN Revision", GetSVNRevision()},
	{"host", GetBuildHostname()},
	{"commandline", util.GetCommandLine() },	
	{"filteringSteps", filteringStep},
	}
	
	flattenTable(LIPRes, "product", stats)
	flattenTable(amgRes, "lastAMG", stats)
	table.append(stats, getAMGStats(amg))
	
	util.printStats(stats)
	if outdir ~= nil  then	
		util.writeFileStats(stats, outdir.."/filtering_stats.txt")
	end

end
	
WriteProfileData("pd.pdxml")


do return end


if true then

v1 = RandomFunction()
smooth(v1, amg, 10, A)
v2 = RandomFunction()
smooth(v2, amg, 10, A)
v3 = RandomFunction()
smooth(v3, amg, 10, A)
VecScaleAssign(v1, 1/VecNorm(v1), v1)
VecScaleAssign(v2, 1/VecNorm(v2), v2)
VecScaleAssign(v3, 1/VecNorm(v3), v3)
--amg:add_watchnode(0,0)
test_ext({v1, v2, v3}, "pevp")
test_ext({v1}, "v1")

end
breakpoint()