--------------------------------------------------------------------------------
--[[!
-- \file apps/amg/famgneumann.lua
--
-- \author Martin Rupp, Andreas Vogel
--
-- \brief Test Script for FAMG with many neumann bnd
--
]]--
--------------------------------------------------------------------------------
-- ProfileLUA(true)
SetOutputProfileStats(false)

ug_load_script("ug_util.lua")
ug_load_script("util.lua")

-- command line parameters
--------------------------------------------------------------------------------------------------

paramOutputProcess = util.GetParamNumber("-outproc", 0, "Output processor")
numRefs = util.GetParamNumber("-numRefs", 1, "number of refinements")
numPreRefs = util.GetParamNumber("-numPreRefs", 0, "number of refinements before ")


maxBase = util.GetParamNumber("-maxBase", 600, "maximum number of nodes before base solver is used")
maxLevels = util.GetParamNumber("-maxLevels", 30, "maximum number of AMG levels")

RAepsilon = util.GetParamNumber("-RAepsilon", 1, "epsilon for rotated anisotropy")
RAalpha = util.GetParamNumber("-RAalpha", 0, "alpha (degrees) for rotated anisotropy")

epsx = util.GetParamNumber("-epsx", 1, "epsilon in X direction")
epsy = util.GetParamNumber("-epsy", 1, "epsilon in Y direction")


bWriteStats = util.HasParamOption("-bWriteStats", "if true, write stats output")
bWriteMat = util.HasParamOption("-writeMatrices", "if true, write matrix output")
bRSAMG = util.HasParamOption("-RSAMG", "if true, use RSAMG, otherwise FAMG") 

bAggressiveCoarsening = util.HasParamOption("-AC", "if true, use aggressive coarsening")
bExternalCoarsening = util.HasParamOption("-XC", "if true, use external coarsening in FAMG")
jacBaseSolver = util.HasParamOption("-jacBaseSolver", "if true, use Linear Jacobi as Base Solver")

famgDelta = util.GetParamNumber("-delta", 0.1, "Interpolation quality F may not be worse than this (F < m_delta)")
famgTheta = util.GetParamNumber("-theta", 0.9, "with multiple parents, discard pairs with m_theta * F > min F.")

minNodes = util.GetParamNumber("-minNodes", 10, "minimum number of nodes on one processor (for agglomeration)")
preferredNodes = util.GetParamNumber("-preferredNodes", 10, "preferred number of nodes on one processor (for agglomeration)")

outdir = util.GetParam("-outdir", "", "output directory")
b3d = util.HasParamOption("-3d", "if true, 3d, otherwise: 2d")
ninePoint = util.HasParamOption("-ninePoint", "for 2d: if true, use ninepoint, otherwise 5pt")
logtofile = util.GetParam("-logtofile", nil, "the file to write logs to")
bWCycle = util.HasParamOption("-Wcycle", "if true, use W-Cycle, otherwise: V-Cycle.")

iTestvectorSmooths = util.GetParamNumber("-iTestvectorSmooths", 2, "nr of testvector smooth for construction of P")
bUseJacobi = util.HasParamOption("-bUseJacobi", "if true, use jacobi as smoother. otherwise use gs")
iSmooths = util.GetParamNumber("-iSmooths", 2, "nr of pre/postsmooths")

--------------------------------------------------------------------------------------------------
util.CheckAndPrintHelp(
"FAMG Script by Martin Rupp"
)
--------------------------------------------------------------------------------------------------


GetLogAssistant():set_output_process(paramOutputProcess)

matrixWritePath=os.getenv("HOME").."/matrices/"
print("Writing matrices to "..matrixWritePath)
if bWriteMat and outdir == "" then
	outdir = matrixWritePath
end

-- constants
if b3d == true then
	dim = 3
else
	dim = 2
end

iSeed = 23

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1));


gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2_four_bnd.ugx")

if logtofile ~= nil then
	GetLogAssistant():enable_file_output(false, "")
	GetLogAssistant():enable_file_output(true, logtofile..GetProcessRank())
	filename = logtofile
else
	if dim == 2 then
		if ninePoint then 
			discName = "9pt"
		else
			discName = "5pt"
		end
	else
		discName = "27pt"
	end
	filename = "uglog_"..discName.."_procs"..GetNumProcesses().."_numRefs"..numRefs
	if bExternalCoarsening then
		filename = filename.."_XC"
	end	
	if bAggressiveCoarsening then
		filename = filename.."_AC"
	end

	GetLogAssistant():enable_file_output(true, outdir..filename..GetProcessRank())
end

print("ugshell "..util.GetCommandLine())
print("SVN Revision: "..GetSVNRevision()..", compiled on "..GetBuildHostname().." on "..GetCompileDate())


print("Parameters: ")
print("    numPreRefs = "..numPreRefs)
print("    numRefs = "..numRefs)
print("    maxBase = "..maxBase)
print("    dim = "..dim)
print("    gridName = "..gridName)
print("    RAepsilon = "..RAepsilon)
print("    RAalpha = "..RAalpha.." degree")
RAalpha = RAalpha * (2*math.pi/360)
print("    RAalpha = "..RAalpha.." grad")
print("    epsx = "..epsx)	
print("    epsy = "..epsy)
print("    gridName = "..gridName)

function writeln(...)
	write(...)
	write("\n")
end

--------------------------------
-- User Data Functions (begin)
--------------------------------
	function cbDirichletBnd2d(x, y, t)
		return true, 0 -- 0 -- 5.9991	
	end
	function cbDirichletBnd3d(x, y, z, t)
		return true, 5.9991		
	end	
	dirchletBnd = LuaCondUserNumber("cbDirichletBnd"..dim.."d")
	
	
	function NeumannNull(x,y,t)
		return 0.0
	end
		
	
	
	function cbSinRhs2d(x, y, t)
		return 0
		--local s = 2*math.pi
		--return s*s*(math.sin(s*x) + math.sin(s*y))
	end
	sinRhs2d = LuaUserNumber("cbSinRhs2d")
	
	
	function cbSinDirichletBnd2d(x, y, t)
		return true, 0
		--local s = 2*math.pi
		--return true, math.sin(s*x) + math.sin(s*y)
	end
	sinDirchletBnd2d = LuaCondUserNumber("cbSinDirichletBnd2d")
	
	-- anisotropic diffusion in corners
	function cbAnisoDiffTensor2d(x, y, t)
		local fac1 = 1/1+exp(alpha*(x+y-1)) 
		local fac2 = 1- fac1
		return	fac1*epsx+fac2, 0, 
				0, fac1+fac2*epsy
	end
--	anisoDiffTensor2d = LuaUserMatrix("cbAnisoDiffTensor2d")

	
	function CreateRotatedAnisotropyMatrix2d(alpha, epsilon)
		local sinalpha = math.sin(alpha)
		local cosalpha = math.cos(alpha)
		RAmat = ConstUserMatrix()
		-- print(sinalpha.." "..cosalpha.." "..epsilon)
		-- print((sinalpha*sinalpha + epsilon*cosalpha*cosalpha)..", "..(1-epsilon)*sinalpha*cosalpha)
		-- print((1-epsilon)*sinalpha*cosalpha..", "..epsilon*sinalpha*sinalpha + cosalpha*cosalpha)
		RAmat:set_entry(0,0,sinalpha*sinalpha + epsilon*cosalpha*cosalpha)
		RAmat:set_entry(0,1,(1-epsilon)*sinalpha*cosalpha)
		RAmat:set_entry(1,0,(1-epsilon)*sinalpha*cosalpha)
		RAmat:set_entry(1,1, epsilon*sinalpha*sinalpha + cosalpha*cosalpha)
		return RAmat		
	end	


	-- hedgehog diffusion 
	function cbHedgehogDiffTensor2d(x, y, t)
	
		if x<0 then 
			if y<0 then return 1.0, 0, 0, 1.0 
			else return epsx, 0, 0, 1.0 end
		else
			if y<0 then return 1.0, 0, 0, epsy 
			else return 1.0, -1.0, -1.0, 1.0 end			
		end
		
		-- should never happen
		return 0, 0, 0, 0 
	end
	
	
	-- discontinuous coefficient
	function cbJumpDiffTensor2d(x, y, t)
		if (math.abs (x)<0.5 and math.abs(y)<0.5) then
			return	epsx, 0, 0, epsy
		else 
			return	1, 0, 0, 1
		end
	end	
	

if dim == 3 then
	diffusionMatrix =  util.CreateConstUserMatrix3d(1, 0, 0,
													0, 1, 0,
													0, 0, 1)
	velocityField = util.CreateConstUserVector3d(0, 0, 0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd

else
	problem = "rotatedAniso"
	-------------------------------------------
	--  Setup User Functions
	-------------------------------------------
	if problem == "rotatedAniso" then
	diffusionMatrix = CreateRotatedAnisotropyMatrix2d(RAalpha, RAepsilon)
	velocityField = ConstUserVector2d(0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd
	end
	
	if problem == "hedgehog" then
	diffusionMatrix = LuaUserMatrix("cbHedgehogDiffTensor2d")
	velocityField = ConstUserVector2d(0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd
	end
	
	if problem == "jump" then
	diffusionMatrix = LuaUserMatrix("cbJumpDiffTensor2d")
	velocityField = ConstUserVector2d(0)
	reaction = ConstUserNumber(0)
	rhs = ConstUserNumber(0)
	dirichlet = dirchletBnd
	end
end

--------------------------------
-- User Data Functions (end)
--------------------------------
tBefore = os.clock()

-- create Instance of a Domain
dom = Domain()
if LoadDomain(dom, gridName) == false then
print("Loading Domain failed.")
exit()
end


-- create Refiner
assert(numPreRefs < numRefs, "numPreRefs must be smaller than numRefs");

refiner = GlobalDomainRefiner(dom)
if numPreRefs > 0 then 
	for i=1,numPreRefs do
	print("refining "..i.."..")
	refiner:refine()
	end
end
-- Distribute the domain to all involved processes
assert(DistributeDomain(dom) == true, "Error while Distributing Grid.")

-- Perform post-refine
if numRefs > 0 then
	for i=numPreRefs+1,numRefs do
	print("refining "..i.."..")
	refiner:refine()
	end
end
tGrid = os.clock()-tBefore

-- get subset handler
sh = dom:subset_handler()


-- write grid to file for test purpose
-- SaveDomain(dom, "refined_grid.ugx")

-- create Approximation Space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

-----------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
-----------------------------------------------------------------
if dim == 2 then
upwind = WeightedUpwind2d(); 
else
upwind = WeightedUpwind3d();
end
upwind:set_weight(0.0)
elemDisc = ConvectionDiffusion("c", "Inner")
elemDisc:set_disc_scheme("fv1")
elemDisc:set_upwind(upwind)
elemDisc:set_diffusion(diffusionMatrix)
elemDisc:set_velocity(velocityField)
elemDisc:set_reaction_rate(reaction)
elemDisc:set_source(rhs)


-----------------------------------------------------------------
--  Setup Neumann Boundary
-----------------------------------------------------------------

if true then
	neumannBND = NeumannBoundary("Inner")
	--neumannBND:add("NeumannNull", "c", "Left,Top,Right")

	-----------------------------------------------------------------
	--  Setup Dirichlet Boundary
	-----------------------------------------------------------------

	
	dirichletBND = DirichletBoundary()
	if util.HasParamOption("-b") then
		dirichletBND:add(dirichlet, "c", "Bottom")
	else
		dirichletBND:add(dirichlet, "c", "Left")
	end
else
	neumannBND = nil
	dirichletBND = DirichletBoundary()
	dirichletBND:add(dirichlet, "c", "Bottom,Top,Left,Right")
end


-------------------------------------------
--  Setup Domain Discretization
-------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
if neumannBND ~= nil then 
	domainDisc:add(neumannBND)
end
domainDisc:add(dirichletBND)


-------------------------------------------
--  Algebra
-------------------------------------------
print ("Setting up Algebra Solver")

-- matrix and vectors
A = MatrixOperator()
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
b2 = GridFunction(approxSpace)

-- set initial value


-- init Operator and set dirichlet values in start iterate
print("Assemble A and b...")
tBefore = os.clock()
domainDisc:assemble_linear(A, b)
tAssemble = os.clock()-tBefore
print ("done")

-- write matrix for test purpose
if bWriteMat then
	SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
	SaveVectorForConnectionViewer(b, "Rhs.vec")
end

-- create algebraic Preconditioners

print ("create preconditioners... ")
jac = Jacobi()
jac:set_damp(0.8)

gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel()
ilu = ILU()
ilut = ILUT()


if bUseJacobi then
	presmoother = jac
	postsmoother = jac
else
	presmoother = sgs
	postsmoother = sgs
end




-- create Base Solver
baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(1000)
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


-- create Convergence Check
convCheck = ConvCheck()
convCheck:set_maximum_steps(100)
convCheck:set_minimum_defect(1e-16)
convCheck:set_reduction(1e-16)

print("done.")
-- create Linear Solver


-- Apply Solver

b2 = GridFunction(approxSpace)
b2:assign(b)

-------------------------------------------
--  Apply Solver
-------------------------------------------
-- 1. init operator


dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_vtk_output(false)

if true then
	-- Gemoetric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(5) -- 1
	gmg:set_base_solver(base)
	gmg:set_smoother(jac)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_debug( dbgWriter)
	transfer = StdTransfer(approxSpace)
	--transfer:set_debug(dbgWriter)
	gmg:set_transfer(transfer)
	
	gmgSolver = LinearSolver()
	gmgSolver:set_preconditioner(gmg)
	
	gmgCGSolver = CG()
	gmgCGSolver:set_preconditioner(gmg)
	
	-- create Convergence Check
	gmgConvCheck = ConvCheck()
	gmgConvCheck:set_maximum_steps(10)
	gmgConvCheck:set_minimum_defect(1e-15)
	gmgConvCheck:set_reduction(1e-20)
	
	gmgSolver:set_convergence_check(gmgConvCheck)
	gmgCGSolver:set_convergence_check(gmgConvCheck)
	
	u = GridFunction(approxSpace)
	b = GridFunction(approxSpace)
	--u:set_random(-1.0, 1.0)
	u:set(4.0)
	
	linOp = AssembledLinearOperator()
	linOp:set_discretization(domainDisc)

	-- 1. init operator
	print("GMG");
	if AssembleLinearOperatorRhsAndSolution(linOp, u, b) == false then 
		print("Could not assemble operator"); exit(); 
	end
	
	srand(iSeed)
	u:set_random(-1.0, 1.0)
	u:set(4.0)	
	domainDisc:adjust_solution(u)
	
		-- 2. apply solver
	print("Apply solver.")
	if ApplyLinearSolver(linOp, u, b, gmgSolver) == false then
		print("Could not apply linear solver.");
	end
	print("done");
	
		-- 1. init operator
	
	SaveVectorForConnectionViewer(u, A, "gmg_u.vec")	
	SaveVectorForConnectionViewer(b, A, "gmg_b.vec")
end
	
