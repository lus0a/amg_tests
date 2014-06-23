-- hello

ug_load_script("laplace_sin_filtering.lua")

checkEVs = util.GetParamNumber("-checkEVs", 8)

if false then
amg = FAMGPreconditioner()
end
--amg:set_delta(0.031)
amg:set_delta(10)
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



--amg:set_testvector_smooths(0)
amg:set_injection_restriction(false)
--------------------------------------------------------------------------------------------

print("\n\n--------------------------------------------")

s = StringTable()
s:set_row_seperator(0, "-")
s:set_row_seperator(1, "-")
s:set_col_seperator(0, "|")
s:set_col_seperator(1, "|")
mytab={}



output = ""
bla = 0
function doTVcheck(myTV)
	local tv = CreateSinVectors(myTV)
	amg:set_testvectors(CreateTV(myTV))
	solver:init(A)
	for i, v in pairs(tv) do
		bla = bla + 1
		SaveVectorForConnectionViewer(v, bla..".vec")
		amg:check_interpolation(v, "__"..bla.."__.vec")
	end
	lis=amg:level_information_str()
	print(lis)
	
	for x=0,checkEVs do
		mytab[x+1] = {}
		s:set(x+1, 0, x)
		s:set(0, x+1, x)
		for y=0,checkEVs do
			avg,last = solve_sin(x,y)			
			mytab[x+1][y+1] = last 
			s:set(x+1, y+1, mytab[x+1][y+1])
		end
	end
	--solve_random()
	s:set_row_seperator(s:num_rows(), "-")
	s:set_col_seperator(s:num_cols(), "|")
	
	print("\n")
	print(util.TableToText(myTV))
	
	output = output.."\n".. "--------------------------------"
	output = output.."\n".. util.TableToText(myTV).."\n"
	print("\n")
	print(s:to_csv("\t"))
	print(s:__tostring())
	print(s:to_latex())
	
	output = output..lis
	output = output.."\n".. s:to_csv("\t") 
	print("\n")
end
if true then
doTVcheck({ {0, 0} })
do return end
doTVcheck({ {1, 1}})

doTVcheck({ {1, 2}})
doTVcheck({ {2, 1}})
doTVcheck({ {2, 2}})
doTVcheck({ {3, 3}})

doTVcheck({ {0, 0}, {1, 1} })
doTVcheck({ {0, 0}, {0, 1} })
doTVcheck({ {0, 0}, {1, 2}} )
doTVcheck({ {1, 1}, {1, 2}} )
doTVcheck({ {1, 2}, {2, 1}} )
doTVcheck({ {2, 1}, {2, 2}} )

doTVcheck({ {0, 0}, {0, 1}, {1, 0} })
doTVcheck({ {0, 0}, {1, 1}, {1, 2} })
doTVcheck({ {1, 1}, {2, 2}, {3, 3} })
doTVcheck({ {1, 1}, {1, 2}, {2, 1} })
end

doTVcheck({ {0, 0} })
doTVcheck({ {1, 1}})
doTVcheck({ {1, 2} })


fout = io.open("filtering_s2_refs_"..numRefs..".txt", "w")
fout:write(util.GetCommandLine())
fout:write(output)
do return end
