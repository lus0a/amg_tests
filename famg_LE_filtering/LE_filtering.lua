



--PrintBuildConfiguration()

ug_load_script("../setup.lua")


big = false
--big = true

if big == true then
	gridName = util.GetParam("-grid", "amgtests/unstructured_33k.ugx")
else
	gridName = util.GetParam("-grid", "amgtests/unstructured_3k.ugx")
end


util.CheckAndPrintHelp("");

matrixWritePath=CreateOutputDirectory("RAMmatrices")

dom = MyLoadGrid({ gridName = gridName, requiredSubsets ={"Inner", "Bottom" }} )

-- needed for unit_square_quads_2x2_4bnd, 
TranslateDomain(dom, 1, 1, 0)
ScaleDomain(dom, 0.5, 0.5, 0.5)
--d = 0.001
--RandomizeDomain(dom, d, d, d)
approxSpace, domainDisc = CreateLEDisc({Domain = dom, BoundaryName = "Bottom"})

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

convCheck = ConvCheck(30, 1e-11, 1e-12, false)
convCheck:set_supress_unsuccessful(true)
precond = Jacobi()




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

--amg:set_matrix_write_path(matrixWritePath)
--amg:set_write_testvectors(true)
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

function myOut(x1, y1, x2, y2)
	local a = CreateSinGridFunction2(x1, y1, x2, y2)
	SaveVectorForConnectionViewer(a, x1..y1..x2..y2..".vec")
end

checkEVs = util.GetParamNumber("-checkEVs", 7)

if false then
amg = FAMGPreconditioner()
end
--amg:set_delta(0.031)
amg:set_delta(0.2)

amg:set_max_nodes_for_base(1)
amg:set_max_levels(2)

numSmooths = 4 --util.GetParamNumber("-amgSmooths", 2)

amg:set_num_presmooth(numSmooths)
amg:set_num_postsmooth(numSmooths)

smoother =  SymmetricGaussSeidel()
--smoother = Jacobi(0.6)
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


matrixWritePath=CreateOutputDirectory("RAMmatrices")
amg:set_matrix_write_path(matrixWritePath)
amg:set_write_testvectors(true)
amg:set_cycle_type(2)
amg:set_Y_cycle(2, 1e-12, 0.2)
--amg:set_testvector_smooths(0)
amg:set_injection_restriction(false)
--amg:add_watchnode(0, 0)

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

test2(TestvectorsLinear2d())
test2(TestvectorsLinearElasticity2d())
test2(TestvectorsConstant())

do return end

amg:set_testvectors(TestvectorsConstant())
for i=1,10 do
	v1 = RandomFunction()
	v2 = RandomFunction()
	v3 = RandomFunction()
	
	
	smooth(v1, amg, 4, A)
	amg:set_one_init(true)
	smooth(v2, amg, 4, A)
	smooth(v3, amg, 4, A)
	
	
	if test_ext({ v1, v2, v3 }, "v123") < 0.3 then break end
	amg:set_one_init(false)
end



do return end
test_ext({ CreateSinGridFunction2(1, 1, 1, 1), CreateSinGridFunction2(0, 1, 0, 1), CreateSinGridFunction2(1, 0, 1, 0) }, "1111, 1010, 0101")
test_ext({ CreateSinGridFunction2(0, 0, 0, 0), CreateSinGridFunction2(0, 1, 0, 1), CreateSinGridFunction2(1, 0, 1, 0) }, "0000, 1010, 0101")
test2(TestvectorsLinear2d())
test2(TestvectorsLinearElasticity2d())
test2(TestvectorsConstant())

test_ext({ CreateSinGridFunction2(1, 1, 1, 1), CreateSinGridFunction2(2, 2, 2, 2), CreateSinGridFunction2(3, 3, 3, 3) }, "1111, 2222, 3333")


test_ext({ CreateSinGridFunction2(0, 0, 0, 0), CreateSinGridFunction2(0, 1, 0, 1), CreateSinGridFunction2(1, 0, 1, 0) }, "0000, 1010, 0101")

test_ext({ CreateSinGridFunction2(1, 1, 1, 1), CreateSinGridFunction2(2, 1, 2, 1), CreateSinGridFunction2(1, 2, 1, 2) }, "1111, 2121, 1212")

test_ext({ CreateSinGridFunction2(1, 1, 1, 1), CreateSinGridFunction2(0, 0, 0, 0) }, "1111, 0000")
test_ext({ CreateSinGridFunction2(1, 1, 1, 1), CreateSinGridFunction2(0, 1, 0, 1) }, "1111, 0101")
--breakpoint()


