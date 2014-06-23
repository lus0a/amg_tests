
--
ug_load_script("laplace_sin_filtering.lua")

amg:set_delta(0.1)
--amg:set_delta(0.1)
amg:set_max_nodes_for_base(1)
amg:set_max_levels(2)

numSmooths = 4 --util.GetParamNumber("-amgSmooths", 2)

amg:set_num_presmooth(numSmooths)
amg:set_num_postsmooth(numSmooths)

smoother = Jacobi(0.6)
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


result3 = {}
result3_2 = {}
scalingTable = {1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6}
function check3(x, fx, fy, scal)
	print("fx = "..fx..", fy = "..fy..", scal = "..scal)
	ev0 = CreateSinGridFunction(0, 0)
	ev2 = CreateSinGridFunction(fx, fy)
	t0 = Eval(ev0 + scal*ev2)

	amg:set_testvectors( TestvectorsExternal({t0}) )
	solver:init(A)
	
	avg,last = solve_sin(0,0)
	result3[x][scal] = last
	avg,last = solve_random()
	result3_2[x][scal] = last	
end	

for x=1,15 do 
	result3[x] = {}
	result3_2[x] = {}
	for i, y in pairs(scalingTable) do
		check3(x, x, x, y)
	end
end	


--util.PrintTable(result3)
print(util.StringTableFromTable(result3, scalingTable):to_csv("\t"))
print(util.StringTableFromTable(result3_2, scalingTable):to_csv("\t"))

