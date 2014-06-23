-- hello

ug_load_script("laplace_sin_filtering.lua")

checkEVs = util.GetParamNumber("-checkEVs", 7)

if false then
amg = FAMGPreconditioner()
end
--amg:set_delta(0.031)
amg:set_delta(0.1)

amg:set_max_nodes_for_base(1)
amg:set_max_levels(2)

numSmooths = 2 --util.GetParamNumber("-amgSmooths", 2)

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

--SetDebugLevel(debugID.
if false then
amg = FAMGPreconditioner()
end


--amg:set_testvector_smooths(0)
amg:set_injection_restriction(false)
--amg:add_watchnode(0, 0)


s = StringTable()
s:set_row_seperator(0, "-")
s:set_row_seperator(1, "-")
s:set_col_seperator(0, "|")
s:set_col_seperator(1, "|")
mytab={}

	
data = ""
for x=0,checkEVs do	
	mytab[x+1] = {}
	s:set(x+1, 0, x)
	s:set(0, x+1, x)
	for y=0,checkEVs do
		
		print("x = "..x.." y = "..y)
		amg:set_testvectors( TestvectorsExternal({ CreateSinGridFunction(x, y) }) )
		--amg:set_testvectors( TestvectorsConstant() )
		solver:init(A)
		avg,last = solve_random()
		mytab[x+1][y+1] = last 
		s:set(x+1, y+1, last)
		data = data..x.." "..y.." "..last.."\n"
	end
end

print("numRefs "..numRefs)
print(s)
print(s:to_csv("\t"))

fout = io.open("output_sm"..numSmooths..".txt", "a")
fout:write("\n\n\n\n\n\n\n")
fout:write(amg:config_string().."\n")
fout:write("numRefs "..numRefs.."\n")
fout:write(s:to_csv("\t").."\n")

do return end
str = [[
set dgrid3d 50,50
set hidden3d
set cntrparam levels 10
set samples 200
set isosamples 210
set xlabel "X axis"
set ylabel "Y axis"
set zlabel "Z axis" offset 1, 0
set view 60, 30, 0.85, 1.1
set title "3D gnuplot demo - contour plot"
set contour
set pm3d interpolate 2,2
]].."splot \"-\" with pm3d \n"
..data
.."e\n"

PipeOpen("gnuplot -persist"):write(str)
breakpoint()






