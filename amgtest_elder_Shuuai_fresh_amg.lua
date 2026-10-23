----------------------------------------------------------
--
--   Lua - Script to test RSAMG for Elder-Problem
--
----------------------------------------------------------
SetOutputProfileStats(false)

ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

-- constants
if util.HasParamOption("-3d") then
	dim = 3
else
	dim = 2
end

gridName	= "../Examples/grids/elder_quads_8x2_bnd.ugx"
numRefs 	= util.GetParamNumber("-numRefs", 6)
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0)
--endTime		= 10 * 365 * 24 * 60 * 60 -- [s] = 10 years
endTime		= 10 * 24 * 60 * 60 -- [s] = 10 days
dt			= 10 * 24 * 60 * 60 -- [s] = 10 days

maxBase = util.GetParamNumber("-maxBase", 20000)
maxLevels = util.GetParamNumber("-maxLevels", 4)
-- Tuned for this coupled Elder problem; each setting can be overridden.
numSmooth = util.GetParamNumber("-numSmooth", 2)
epsilonStrong = util.GetParamNumber("-epsilonStrong", 0.05)
fSmoothing = util.GetParamBool("-fsmoothing", false)
matrixPath = util.GetParam("-matrixPath", "") -- Nonempty enables matrix dumps.
smootherName = util.GetParam("-smoother", "ilu")
precondName = util.GetParam("-precond", "amg")

vtk_file_name = "amgtestvtk/amgtest_elder_optimized"

bFileOutput = false
bOutput = true

bRSAMG = true


Porosity = 0.1
MolecularDiffusion = 3.565e-6
Permeability = 4.845e-13
Viscosity = 1e-3

function PressureStart(x, y, t, si)
	return 9810 * (150 - y)
end

function ConcentrationDirichletBnd(x, y, t)
		if x > 150 and x < 450 then
			return true, 1.0
		end

	    return false, 0.0
end

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 2));

-- create the domain, load the grid and refine it
dom = util.CreateDomain(gridName, numPreRefs, {"Inner", "LeftRightBnd", "BottomBnd", "TopBnd","TopCorners"})
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:add_fct("p", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()

-- Order the DoFs:
OrderLex (approxSpace, "y")

--------------------------------------------------------------------------------
--  Setup the spatial FV Discretization
--------------------------------------------------------------------------------

Gravity = ConstUserVector(0.0)
Gravity:set_entry(1, -9.81)

-- create dirichlet boundary for concentration
dirichletBND = DirichletBoundary()
dirichletBND:add("ConcentrationDirichletBnd", "c", "TopBnd")
dirichletBND:add(0, "c", "BottomBnd")
-- TopCorners
dirichletBND:add(0, "p", "TopCorners")

-- density
function DensityFct(c) return 1000 + 200 * c end
function DDensityFct_c(c) return 200 end
Density = LuaUserFunctionNumber("DensityFct", 1);
Density:set_deriv(0, "DDensityFct_c");

-- Darcy Velocity (as a linker)
DarcyVelocity = DarcyVelocityLinker();
DarcyVelocity:set_permeability(Permeability)
DarcyVelocity:set_viscosity(Viscosity)
DarcyVelocity:set_density(Density)
DarcyVelocity:set_gravity(Gravity)

-- flow equation
FlowEq = ConvectionDiffusionFV1("p", "Inner")
FlowEq:set_mass_scale(0.0)
FlowEq:set_mass(Porosity * Density)
FlowEq:set_flux(DarcyVelocity * Density)

-- transport equation
TransportEq = ConvectionDiffusionFV1("c", "Inner")
TransportEq:set_upwind(FullUpwind())
TransportEq:set_mass_scale(Porosity * Density)
TransportEq:set_velocity(DarcyVelocity * Density)

-- molecular Diffusion
Diffusion = ScaleAddLinkerMatrix()
Diffusion:add(Porosity, MolecularDiffusion)
TransportEq:set_diffusion(Diffusion* Density)

-- associate the density with the unknown of the transport equation (i.e. c)
Density:set_input(0, TransportEq:value())
-- associate the Darcy velocity with the unknown of the flow equation (i.e. p)
DarcyVelocity:set_pressure_gradient(FlowEq:gradient())

-- set up the global discretization
domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(TransportEq)
domainDisc:add(FlowEq)
domainDisc:add(dirichletBND)

-----------------------------------------------------------------------------------------
-- Ploter for the results
------------------------------------------------------------------------------------------

out = VTKOutput()
out:clear_selection()
out:select_nodal("c", "c")
out:select_nodal("p", "p")
out:select(DarcyVelocity, "q")

-- grid function for the solution
u = GridFunction(approxSpace)
--b = GridFunction(approxSpace)
--b2 = GridFunction(approxSpace)

-- Set the initial condition
Interpolate("PressureStart", u, "p")
Interpolate(0, u, "c")

--------------------------------------------------------------------------------
--  Set up the solver
--------------------------------------------------------------------------------

--[[
-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)

-- init Operator
linOp:init_op_and_rhs(b)
linOp:set_dirichlet_values(u)
--]]

----[[
jac = Jacobi()
jac:set_damp(0.6)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()

if smootherName == "gs" then
	smoother = gs
elseif smootherName == "sgs" then
	smoother = sgs
elseif smootherName == "ilu" then
	smoother = ilu
elseif smootherName == "ilut" then
	smoother = ilut
else
	error("Unsupported -smoother '"..smootherName.."' (use gs, sgs, ilu, or ilut)")
end


-- create Base Solver
baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-16)
baseConvCheck:set_reduction(1e-16)
baseConvCheck:set_verbose(false)


base = LU()

print ("create AMG... ")
amg = RSAMGPreconditioner()
--amg:enable_aggressive_coarsening_A(2)

vectorWriter = GridFunctionPositionProvider()
vectorWriter:set_reference_grid_function(u)
amg:set_position_provider(vectorWriter)

if matrixPath ~= "" then
	amg:set_matrix_write_path(matrixPath)
end

amg:set_num_presmooth(numSmooth)
amg:set_num_postsmooth(numSmooth)
amg:set_cycle_type(1)
amg:set_presmoother(smoother)
amg:set_postsmoother(smoother)

amg:set_base_solver(base)

amg:set_max_levels(maxLevels)
--amg:set_max_levels(numRefs)

-- amg:set_min_nodes_on_one_processor(200) not functional
amg:set_max_nodes_for_base(maxBase)
amg:set_max_fill_before_base(0.7)
amg:set_fsmoothing(fSmoothing)
amg:set_epsilon_strong(epsilonStrong)
--amg:set_epsilon_truncation(0)

--amg:set_one_init(true)  -- EXPERIMENT: fresh AMG setup on every Newton step

print("create AMG done.")
----]]



--[[
-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(30)
--convCheck:set_maximum_steps(300)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)
--convCheck:set_reduction(1e-2)

-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(amg)
linSolver:set_convergence_check(convCheck)
print("create linSolver done.")

log = GetLogAssistant();

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init
linOp:init_op_and_rhs(b2)
linSolver:init(linOp)

print(amg:level_information_str())

-- 2. apply solver
tBefore = os.clock()
linSolver:apply_return_defect(u,b)
tSolve = os.clock()-tBefore
WriteGridFunctionToVTK(u, "Solution")


print("Solver done")

if bOutput then
WriteGridFunctionToVTK(u, "Solution")
end

--amg:check(u, b)  --not work
amg:check2(u, b)

log:set_debug_levels(0)

printf = function(s,...)
	print(s:format(...))
end -- function
----]]

--------------------------------------------------------------------------------
--  Set up the solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace
if precondName == "amg" then
	precondDesc = {type = "amg", amg = amg}
elseif precondName == "ilu" or precondName == "ilut" then
	precondDesc = {type = precondName}
else
	error("Unsupported -precond '"..precondName.."' (use amg, ilu, or ilut)")
end

solverDesc = {
	type = "newton",

	convCheck = {
		type		= "standard",
		iterations 	= 30,		-- maximum number of iterations
		absolute	= 5e-6,		-- absolut value of defect to be reached; usually 1e-7 - 1e-9
		reduction	= 1e-10,	-- reduction factor of defect to be reached; usually 1e-6 - 1e-8
		verbose		= true			-- print convergence rates if true
	},

	linSolver =
	{
		type = "bicgstab",
		--[[
		precond = {
	        type            = "gmg",        -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
	        smoother		= "ilu", --"ilu",		-- pre- and postsmoother, only for gmg ["ilu", "ilut", "jac", "gs", "sgs"]
	        cycle           = "V",          -- gmg-cycle ["V", "F", "W"]
	        preSmooth       = 2,            -- number presmoothing steps
	        postSmooth      = 2,            -- number postsmoothing steps
	        rap             = false,        -- comutes RAP-product instead of assembling if true
	        baseLevel       = numPreRefs,	-- gmg - baselevel
	        baseSolver		= "lu",
		},
		--]]
		precond = precondDesc,

		convCheck = {
	        type            = "standard",
	        iterations      = 500,          -- number of iterations
	        absolute        = 1e-8,        -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be larger than in newton section)
	        reduction       = 1e-5,         -- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be larger than in newton section)
	        verbose         = true,         -- print convergence rates if true
		}
	}
}

solver = util.solver.CreateSolver(solverDesc)


------------------------------------------------------------------------------------------
-- Set up and apply the time stepping scheme
------------------------------------------------------------------------------------------

-- Set up the time discretization and the time stepping scheme
timeDisc = ThetaTimeStep(domainDisc, 1.0) -- Implicit Euler: 1.0
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name, out))

-- Compute the time steps
timeIntegrator:set_time_step(dt)
local solveStart = os.clock()
assert(timeIntegrator:apply(u, endTime, u, 0.0), "Elder time integration failed")
print(string.format("Time integration CPU: %.3f s", os.clock() - solveStart))

------------------------------------------------------------------------------------------
-- End of File
------------------------------------------------------------------------------------------

--[[
-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- create operator from discretization
linOp = AssembledLinearOperator()
linOp:set_discretization(domainDisc)

-- get grid function
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
b2 = GridFunction(approxSpace)

-- set initial value
u:set_random(-1.0, 1.0)

-- init Operator
linOp:init_op_and_rhs(b)
linOp:set_dirichlet_values(u)

jac = Jacobi()
jac:set_damp(0.6)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()


-- create Base Solver
baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(500)
baseConvCheck:set_minimum_defect(1e-16)
baseConvCheck:set_reduction(1e-16)
baseConvCheck:set_verbose(false)

if false then
	base = LinearSolver()
	base:set_convergence_check(baseConvCheck)
	base:set_preconditioner(jac)
else
	base = LU()
end

-- create AMG ---
-----------------

print ("create AMG... ")
amg = RSAMGPreconditioner()
--amg:enable_aggressive_coarsening_A(2)


vectorWriter = GridFunctionPositionProvider()
vectorWriter:set_reference_grid_function(u)
amg:set_position_provider(vectorWriter)


if bOutput then
amg:set_matrix_write_path("/home/lus0a/amgMatrices/")
end

amg:set_num_presmooth(3)
amg:set_num_postsmooth(3)
amg:set_cycle_type(1)
--amg:set_presmoother(gs)
--amg:set_postsmoother(gs)
amg:set_presmoother(ilu)
amg:set_postsmoother(ilu)

amg:set_base_solver(base)

amg:set_max_levels(4)

-- amg:set_min_nodes_on_one_processor(200) not functional
amg:set_max_nodes_for_base(maxBase)
amg:set_max_fill_before_base(0.7)
amg:set_fsmoothing(true)
--amg:set_epsilon_truncation(0)

--amg:set_one_init(true)  -- EXPERIMENT: fresh AMG setup on every Newton step

--amg:tostring()

-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(30)
convCheck:set_minimum_defect(1e-11)
convCheck:set_reduction(1e-12)

print("done.")
-- create Linear Solver
linSolver = LinearSolver()
linSolver:set_preconditioner(amg)
linSolver:set_convergence_check(convCheck)

-- Apply Solver

log = GetLogAssistant();

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init
linOp:init_op_and_rhs(b2)
linSolver:init(linOp)

print(amg:level_information_str())

-- 2. apply solver
tBefore = os.clock()
linSolver:apply_return_defect(u,b)
tSolve = os.clock()-tBefore
WriteGridFunctionToVTK(u, "Solution")


print("done")

if bOutput then
WriteGridFunctionToVTK(u, "Solution")
end

--amg:check(u, b)  --not work
amg:check2(u, b)

log:set_debug_levels(0)

printf = function(s,...)
	print(s:format(...))
end -- function


if bFileOutput then
output = io.open("output.txt", "a")

output:write("dim\t"..dim)
output:write("\tprocs\t"..NumProcs())
-- output:write("\tgridName\t"..gridName)
-- output:write("\tnumPreRefs\t"..numPreRefs)
output:write("\tnumRefs\t".. numRefs)
output:write("\tanisotropy\t".. RAepsilon)
output:write("\talpha\t"..(RAalpha * (180/math.pi)))
output:write("\tsteps\t"..convCheck:step())
output:write("\tlast reduction\t".. convCheck:defect()/convCheck:previous_defect())
output:write("\tamg_timing_whole_setup\t".. amg:get_timing_whole_setup_ms())
output:write("\tc_A\t".. amg:get_operator_complexity())
output:write("\tc_G\t".. amg:get_grid_complexity())
output:write("\tusedLevels\t"..amg:get_used_levels())
output:write("\ttSolve [ms]\t".. tSolve)
output:write("\n")
print(s)
else

print("AMG Information:")
printf("Setup took %.2f ms",amg:get_timing_whole_setup_ms())
printf("Coarse solver setup took %.2f ms",amg:get_timing_coarse_solver_setup_ms())
printf("Operator complexity c_A = %.2f", amg:get_operator_complexity())
printf("Grid complexity c_G = %.2f", amg:get_grid_complexity())

LI = amg:get_level_information(0)
printf("Level 0: number of nodes: %d. fill in %.2f%%, nr of interface elements: %d (%.2f%%)", LI:get_nr_of_nodes(), LI:get_fill_in()*100, LI:get_nr_of_interface_elements(), LI:get_nr_of_interface_elements()/LI:get_nr_of_nodes()*100)
for i = 1, amg:get_used_levels()-1 do
	LI = amg:get_level_information(i)
	printf("Level %d: creation time: %.2f ms. number of nodes: %d. fill in %.2f%%. coarsening rate: %.2f%%, nr of interface elements: %d (%.2f%%)", i,
		LI:get_creation_time_ms(),  LI:get_nr_of_nodes(), LI:get_fill_in()*100,
		LI:get_nr_of_nodes()*100/amg:get_level_information(i-1):get_nr_of_nodes(), LI:get_nr_of_interface_elements(), LI:get_nr_of_interface_elements()/LI:get_nr_of_nodes()*100)

end


if GetProfilerAvailable() == true then

	create_levelPN = GetProfileNode("c_create_AMG_level")

	function PrintParallelProfileNode(name)
		pn = GetProfileNode(name)
		t = pn:get_avg_total_time_ms()/to100 * 100
		tmin = ParallelMin(t)
		tmax = ParallelMax(t)
		printf("%s:\n%.2f %%, min: %.2f %%, max: %.2f %%", name, t, tmin, tmax)
	end

	if create_levelPN:is_valid() then
		if true then
			print(create_levelPN:call_tree())
			print(create_levelPN:child_self_time_sorted())
		end
		to100 = create_levelPN:get_avg_total_time_ms()
		PrintParallelProfileNode("create_OL2_matrix")
		PrintParallelProfileNode("CalculateTestvector")
		PrintParallelProfileNode("CreateSymmConnectivityGraph")
		PrintParallelProfileNode("calculate_all_possible_parent_pairs")
		PrintParallelProfileNode("color_process_graph")
		PrintParallelProfileNode("FAMG_recv_coarsening_communicate")
		PrintParallelProfileNode("update_rating")
		PrintParallelProfileNode("precalculate_coarsening")
		PrintParallelProfileNode("send_coarsening_data_to_processes_with_higher_color")
		PrintParallelProfileNode("communicate_prolongation")
		PrintParallelProfileNode("create_new_index")
		PrintParallelProfileNode("create_parent_index")
		PrintParallelProfileNode("create_interfaces")
		PrintParallelProfileNode("create_fine_marks")
		PrintParallelProfileNode("FAMGCreateAsMultiplyOf")
		PrintParallelProfileNode("CalculateNextTestvector")
	end

else
	print("Profiler not available.")
end


end
-- print("-----------------------------------------------")
-- amg:check2(u, b);

--]]
