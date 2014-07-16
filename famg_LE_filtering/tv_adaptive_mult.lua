



--PrintBuildConfiguration()

ug_load_script("../setup.lua")


big = false
--big = true

if util.HasParamOption("-33k") then
	gridName = util.GetParam("-grid", "amgtests/unstructured_33k.ugx")
elseif util.HasParamOption("-3k") then
	gridName = util.GetParam("-grid", "amgtests/unstructured_3k.ugx")
elseif util.HasParamOption("-special") then
	gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_4bnd2.ugx")
else
	gridName = util.GetParam("-grid", "unit_square/unit_square_quads_2x2_4bnd.ugx")
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
function createAMG1()
	local amg = FAMGPreconditioner()
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
	amg:set_write_testvectors(true)
	
	return amg
end

function createLS(amg)
	-- create Linear Solver
	linSolver = LinearSolver(amg, convCheck)
	-- choose some solver
	solver = linSolver
end
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

function createAMG()
	amg = createAMG1()
	if false then
	amg = FAMGPreconditioner()
	end
	--amg:set_delta(0.031)
	amg:set_delta(0.3)
	amg:set_max_nodes_for_base(1000)
	amg:set_max_levels(2)
	
	numSmooths = 4 --util.GetParamNumber("-amgSmooths", 2)
	
	amg:set_num_presmooth(numSmooths)
	amg:set_num_postsmooth(numSmooths)
	
	--smoother = Jacobi(0.6)
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




convCheck:set_verbose(true)


lip = LinearIteratorProduct()
lsLIP = LinearSolver(lip, convCheck)

numVecs = util.GetParamNumber("-numVecs", 2)

v={}
for k=1,numVecs do
	v[k]= RandomFunction()
	smooth(v[k], SymmetricGaussSeidel(), 10, A)
end


for filteringStep=1,20 do
	print("--------------------------------------------------------")
	print("                        "..filteringStep)
	print("--------------------------------------------------------")
	
	

	amg = createAMG()
	amg:set_one_init(true)
	createLS(amg)
	
	if util.HasParamOption("-writeMatrices") then
		for k=1,numVecs do
			SaveVectorForConnectionViewer(v[k], "v_"..filteringStep.."."..k..".vec")
		end
		amg:set_matrix_write_path(CreateOutputDirectory("RAMmatrices/"..filteringStep))
	end
	
	amg:set_testvectors(TestvectorsExternal(v))
	solver:init(A)

	lip:add_iterator(amg)
	--aaa=FAMGPreconditioner()
	print(amg:level_information_str())
	
	lsLIP:init(A)
	
	if filteringStep ~= 1 then amgRes = solve_random_d(solver, approxSpace) end
	LIPRes = solve_random_d(lsLIP, approxSpace)	
	if LIPRes.lastReduction < 0.1 then break end
	
	for k=1,numVecs do
		v[k] = RandomFunction()
		smooth(v[k], lip, 10, A)
		VecScaleAssign(v[k], 1/VecNorm(v[k]), v[k])
	end
	
	for k1=1,numVecs do
		for k2=k1+1,numVecs do
			Assign(v[k1], v[k1]-VecProd(v[k1], v[k2])*v[k2])
		end
		VecScaleAssign(v[k1], 1/VecNorm(v[k1]), v[k1])
	end
	
	
	for k1=1,numVecs do
	for k2=1,numVecs do
	print(VecProd(v[k1], v[k2]))
	end end
	
end

if ProcRank() == 0 then
	stats = {
	
	{ "procs", NumProcs() },
	{ "Type", "LE" },
	{ "numVecs", numVecs},
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

