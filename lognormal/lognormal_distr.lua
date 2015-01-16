
--PrintBuildConfiguration()

ug_load_script("ug_util.lua")
ug_load_script("solver_util/setup_rsamg.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("../setup.lua")





--------------------------------------------------------------------------------
-- Checking for parameters (begin)
--------------------------------------------------------------------------------
-- Several definitions which can be changed by command line parameters
-- space dimension and grid file:
dim = util.GetParamNumber("-dim", 2, "Dimension of the problem", {2,3})

if dim == 2 then
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_tri_2x2.ugx")
	gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_8x8.ugx")
	--gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
	zWriteOffset = 0.1
end
if dim == 3 then
	gridName = util.GetParam("-grid", "unit_square/unit_cube_hex.ugx")
	--gridName = "unit_square/unit_cube_tets_regular.ugx"
	zWriteOffset = 3
end

-- refinements:
numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before parallel distribution")
numRefs    = util.GetParamNumber("-numRefs",    3, "number of refinements")

ug_assert(numPreRefs <= numRefs, "It must be choosen: numPreRefs <= numRefs")

-- if withRedist == true, you may alter the redistribution behavior over the
-- parameters declared in util/load_balancing_util.lua.
-- Note that numPreRefs is ignored in this case.
withRedist = util.HasParamOption("-withRedist", false, "Enables a more sophisticated distribution approach.")
staticRedist = util.HasParamOption("-staticRedist", false, "Static distribution approach.")

saveGrids = util.HasParamOption("-saveGrids", false)

--! amount of output
--! set to 0 i.e. for time measurements,
--! >= 1 for writing matrix files etc.
bVerb = util.HasParamOption("-verb", "verbosity. Set for output when no time measurements")		

activateDbgWriter = 0	  
--! set to 0 i.e. for time measurements,
--! >= 1 for debug output: call 'set_debug(dbgWriter)'
--! for the main solver ('gmg')
activateDbgWriter = util.GetParamNumber("-dbgw", 0, "debug writer: 0 for no output.") 

-- controls for gmg						    
numPreSmooth = util.GetParamNumber("-numPreSmooth", 3)
numPostSmooth = util.GetParamNumber("-numPostSmooth", 3)
if util.HasParamOption("-numSmooth") then
numPreSmooth = util.GetParamNumber("-numSmooth", 3)
numPostSmooth = numPreSmooth
end
smooth = util.GetParam("-smooth", "jac", "smoother", {"jac","ilu","gs","sgs"})
baseLev = util.GetParamNumber("-baseLev", 0)
cycleType =  util.GetParam("-cycle", "V", "gmg-cycle type", {"V","W","F"})
bRAP = util.HasParamOption("-rap", "use rap product as level matrices")
bSmoothRim = util.HasParamOption("-smoothRim", "Smooth on Surface rim in gmg")
bEmulateFullGrid = util.HasParamOption("-emulateFullGrid", "use full-refined grid emulation")
bOverlap = util.HasParamOption("-overlap", "overlap communication and computation")

disc = util.GetParam("-disc", "fv1", "Disc type", {"fv1","fv","fe","fecr","fvcr"})
p = util.GetParamNumber("-p", 1, "Order of space")

util.CheckAndPrintHelp("LognormalDist");

--------------------------------------------------------------------------------
-- Checking for parameters (end)
--------------------------------------------------------------------------------

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

	requiredSubsets = {"Inner", "Boundary"}

if withRedist == true then
	dom = util.CreateDomain(gridName, 0, requiredSubsets)
	balancer.RefineAndRebalanceDomain(dom, numRefs)
elseif staticRedist == true then
	balancer.ParseParameters()
	dom = util.CreateDomain(gridName, 0, requiredSubsets)
	balancer.RefineAndRebalanceDomain(dom, numRefs)
else
	dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, requiredSubsets)
end

function Source2d(x, y, t)
	local s = 2*math.pi
	return s*s*(math.sin(s*x) + math.sin(s*y))
end

function DirichletValue2d(x, y, t)
	local s = 2*math.pi
	return true, 0 --math.sin(s*x) + math.sin(s*y)
end

function Source3d(x, y, z, t)
	return 0;
end

function DirichletValue3d(x, y, z, t)
	return true, x
end

	
function do_stuff(rm_nodes, rm_mean_f, rm_sigma_f, rm_sigma)

srand(4)
	--srand(GetClockS()*100)
	Normal_Field=LognormalRandomField(rm_nodes, rm_mean_f, rm_sigma_f, rm_sigma)
	print(Normal_Field:config_string())
	--Normal_Field:set_no_exp()

	
	-- create Approximation Space
	print("Create ApproximationSpace")
	approxSpace = ApproximationSpace(dom)
	approxSpace:add_fct("c", "Lagrange", p)
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	
	--OrderCuthillMcKee(approxSpace, true);
	
	--------------------------------------------------------------------------------
	--  Assembling
	--------------------------------------------------------------------------------
	elemDisc = ConvectionDiffusion("c", "Inner", disc)
	--elemDisc:set_diffusion(1)
	elemDisc:set_diffusion(Normal_Field)
	elemDisc:set_source("Source"..dim.."d")
	
	dirichletBND = DirichletBoundary()
	dirichletBND:add("DirichletValue"..dim.."d", "c", "Boundary")
	
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(elemDisc)
	domainDisc:add(dirichletBND)
	
	--------------------------------------------------------------------------------
	--  Algebra
	--------------------------------------------------------------------------------
	print ("Setting up Algebra Solver")
	
	
	amg = FAMGPreconditioner()
	amgsetup.SetFAMGCommandLineParameters(amg, {
		smooths=2,
		maxNodesForBase=2000,
		delta=0.1	
		})
		
	
	if util.HasParamOption("-writeMatrices") then
		matrixWritePath=CreateOutputDirectory("RAMmatrices")
		amg:set_matrix_write_path(matrixWritePath)
		amg:set_write_f_values(true)
		amg:set_write_testvectors(true)
	end
	
		
	
	-- Geometric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(baseLev)
	gmg:set_base_solver(LU())
	gmg:set_gathered_base_solver_if_ambiguous(true)
	if 	    smooth == "ilu" then gmg:set_smoother(ILU());
	elseif 	smooth == "jac" then gmg:set_smoother(Jacobi(0.8));
	elseif 	smooth == "gs"  then gmg:set_smoother(GaussSeidel());
	elseif 	smooth == "sgs" then gmg:set_smoother(SymmetricGaussSeidel());
	else print("Smoother not set, use -smooth option [ilu, jac, gs, sgs]"); exit(); end
	gmg:set_num_presmooth(numPreSmooth)
	gmg:set_num_postsmooth(numPostSmooth)
	gmg:set_smooth_on_surface_rim(bSmoothRim)
	gmg:set_cycle_type(cycleType)
	gmg:set_rap(bRAP)
	gmg:set_emulate_full_refined_grid(bEmulateFullGrid)
	gmg:set_comm_comp_overlap(bOverlap)
	if activateDbgWriter >= 1 then
		dbgWriter = GridFunctionDebugWriter(approxSpace)
		dbgWriter:set_vtk_output(true)
		gmg:set_debug(dbgWriter)
	end
	
	-- choose precond
	if util.HasParamOption("-amg") then
		print(amg:config_string())
		precond = amg
	else
		precond = gmg
	end
		
	-- create Convergence Check
	convCheck = ConvCheck(100, 1e-99, 1e-12, true)
	
	-- create Linear Solver
	solver = LinearSolver(precond, convCheck)
	--------------------------------------------------------------------------------
	-- Solve
	--------------------------------------------------------------------------------
	
	-- 1. Create Algebra
	A = MatrixOperator()
	u = GridFunction(approxSpace)
	b = GridFunction(approxSpace)
	
	-- 2. Reset start solution
	u:set(0.0)
	domainDisc:adjust_solution(u)
	
	-- 3. assemble matrix and rhs
	domainDisc:assemble_linear(A, b)
	
	if  bVerb == true then
		print("Save matrix and rhs for connection viewer")
		matFile = "Stiffness_"..counter..".mat"
		SaveMatrixForConnectionViewer(u, A, matFile)
		os.execute("mv "..matFile.." "..matFile..".vec") 
		SaveVectorForConnectionViewer(b, "Rhs_"..counter..".vec")
	end
	
	-- 4. init solver for linear Operator
	solver:init(A)
	
	-- 5. apply solver
	solver:apply_return_defect(u,b)
	local lastReduction = convCheck:defect()/convCheck:previous_defect()
	print(lastReduction)
	
	if util.HasParamOption("-ev") then
		nev = 10
		eig = EigenSolver()
		eig:set_linear_operator_A(A)
		eig:set_max_iterations(20)
		eig:set_precision(1e-5)
		eig:set_preconditioner(gmg) --ILUT(0.0))
		eig:set_pinvit(2)
		eig:set_additional_eigenvectors_to_keep(8)
		ev = {}
		for i=1,nev do
			print("adding ev "..i)
			ev[i] = GridFunction(approxSpace)
			if ev0 ~= nil then 
				Assign(ev[i], 1.0*ev0[i])
			else
				ev[i]:set_random(-1.0, 1.0)
			end
			domainDisc:adjust_solution(ev[i])
			eig:add_vector(ev[i])
		end
		
		if ev0 == nil then
			ev0 = {}
			for i=1,nev do
				ev0[i] = ev[i]
			end
		end
		
		eig:apply()
		
		for i=1,nev do
			--WriteGridFunctionToVTK(ev[i], "C"..counter..ev_"..i)
			
			SaveVectorForConnectionViewer(ev[i], "C"..counter.."ev_"..i.."_"..eig:get_eigenvalue(i-1)..".vec")
			
		end
		counter = counter+1
	end
end
counter=0
--rm_nodes = util.GetParamNumber("-RF_N", 100)
--rm_sigma = util.GetParamNumber("-RF_sigma", 0.01)
--rm_sigma_1 = util.GetParamNumber("-RF_sigma1", rm_sigma)
--rm_sigma_2 = util.GetParamNumber("-RF_sigma2", rm_sigma)

--do_stuff(5, 0, 1, 0.01)
do_stuff(11, 0, 1, 0.1)

--do_stuff(10, 0, 10, 0.01)
do return end

sigma= 2
for i=1,14 do
do_stuff(200, 0, 4, sigma)
print("---- sigma = "..sigma.." ----")
sigma = sigma / 2
end

-- end group app_convdiff
--[[!
\}
]]--
