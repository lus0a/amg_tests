-------------------------------------------------------------------------------
--[[!
-- \file amg/laplace_sin_filtering.lua
-- \{
-- We set a Laplace-problem up here
--
]]--
--------------------------------------------------------------------------------


--PrintBuildConfiguration()

ug_load_script("setup.lua")

InitSignals()

dim = util.GetParamNumber("-dim", 2)

if dim == 3 then
gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_4x4x4.ugx")
dom = MyLoadGrid({ dim=3, gridName = gridName, requiredSubsets = requiredSubsets} )
else
gridName = util.GetParam("-grid", "unit_square/unit_square_quads_8x8.ugx")
dom = MyLoadGrid({ dim=2, gridName = gridName, requiredSubsets = requiredSubsets} )

end


requiredSubsets = {"Inner", "Boundary"}
boundaries = {"Boundary"}

util.CheckAndPrintHelp("");


approxSpace = CreateStdApproxSpace(dom)
domainDisc = CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries})

	

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
write("Assembling...")
domainDisc:assemble_linear(A, b)
SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
print(" done.")

amg = FAMGPreconditioner()
SetFAMGCommandLineParameters(amg, {
	--delta=0.030473,
	delta=0.1,
	smooths=4,
	maxNodesForBase=2000	
	})
amg:set_deficit_solve_tolerance(1e-15)

sgs = SymmetricGaussSeidel()
amg:set_presmoother(sgs)
amg:set_postsmoother(sgs)

matrixWritePath=CreateOutputDirectory("RAMmatrices")
amg:set_matrix_write_path(matrixWritePath)
amg:set_write_testvectors(true)

amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
amg:set_testvector_smooths(0)
amg:set_testvector_smoother(Jacobi(0.6))


amg:set_write_f_values(true)
--amg:set_debug_write_all_vectors(true)
amg:set_fsmoothing(true)
amg:set_testvectors(TestvectorsLinear2d())
--amg:set_testvectors(TestvectorsConstant())
amg:set_verbose(5)
--amg = RSAMGPreconditioner()

amg = RSAMGPreconditioner()
convCheck = ConvCheck(30, 1e-12, 1e-10)
linSolver = LinearSolver(amg, convCheck)
solver = linSolver

print(solver:config_string())
print("Init Solver...")

--SetDebugLevel(debugID.AMG, 5)
--SetDebugLevel(debugID.AMG.FAMG.Interpolation, 0)

solver:init(A)

lis=amg:level_information_str()
	print(lis)
print("Solve...")


solve_random(solver, approxSpace)

amg:set_one_init(true)

function solve_lin(linSolver, amg, name)
	linSolver:set_convergence_check(convCheck)	
	linSolver:set_preconditioner(amg)
	linSolver:init(A)
	
	solve_random(linSolver, approxSpace)
	
end

function solve_amg(smoother, numSmooths, name)


	amg:set_presmoother(smoother)
	amg:set_postsmoother(smoother)
	amg:set_num_presmooth(numSmooths)
	amg:set_num_postsmooth(numSmooths)
	
	cgSolver = BiCGStab()
	cgSolver:set_preconditioner(amg)	
	amg:init_smoothers()
	
	print(name..": "..numSmooths.."x "..smoother:config_string())
	
	
	--ScaleDomain(dom, MakeVec(0,0,0), MakeVec(1/stretchX, 1/stretchY, 1.0))
	--SaveVectorForConnectionViewer(u, "amg_u_linear.vec")
	--ScaleDomain(dom, MakeVec(0,0,0), MakeVec(stretchX, stretchY, 1.0))
		
	print("CG, "..name)
	solve_lin(CG(), amg, "CG("..name..")")
	print("LinearSolver, "..name)
	solve_lin(LinearSolver(), amg, "LS("..name..")")	
	

end

solve_amg(SymmetricGaussSeidel(), 2, "SGSx02")
solve_amg(SymmetricGaussSeidel(), 32, "SGSx32")
