ug_load_script("laplace_sin_filtering.lua")




function check2(fx, fy, scal)
	print("fx = "..fx..", fy = "..fy..", scal = "..scal)
	ev0 = CreateSinGridFunction(1, 1)
	ev2 = CreateSinGridFunction(fx, fy)
	t0 = Eval(ev0 + scal*ev2)

	amg:set_testvectors( TestvectorsExternal({t0}) )
	solver:init(A)
	for x=1,checkEVs do
		mytab[x] = {}
		s:set(x, 0, x)
		s:set(0, x, x)
		for y=1,checkEVs do
			avg,last = solve_sin(x,y)			
			mytab[x][y] = last 
			s:set(x, y, mytab[x][y])
		end
	end
	print(s:to_csv("\t"))
end	

s = StringTable()
s:set_row_seperator(0, "-")
s:set_row_seperator(1, "-")
s:set_col_seperator(0, "|")
s:set_col_seperator(1, "|")
mytab={}

check2(4, 4, 0.000001)
check2(4, 4, 0.00001)
check2(4, 4, 0.0001)
check2(4, 4, 0.001)
check2(4, 4, 0.01)
check2(4, 4, 0.1)
check2(4, 4, 0.2)
check2(4, 4, 0.4)
check2(4, 4, 0.8)

do return end
check2(5, 5, 0.1)
check2(5, 5, 0.2)
check2(5, 5, 0.4)
check2(5, 5, 0.8)
check2(5, 5, 1.0)
check2(5, 5, 1.1)
check2(5, 5, 1.2)
check2(5, 5, 1.4)
