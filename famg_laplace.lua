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
outfile = util.GetParam("-outfile", nil, "output directory")


if true then

	requiredSubsets = {"Inner", "Boundary"}
	
	if util.HasParamOption("-LE") then		
		cpu = dim
	else
		cpu = 1
	end
	
	if dim == 3 then
		gridName = util.GetParam("-grid", "unit_square_01/unit_cube_01_hex_4x4x4.ugx")
		dom = MyLoadGrid({ dim=3, gridName = gridName, requiredSubsets = requiredSubsets} )
		boundaries = {"Boundary"}
	else
		if util.HasParamOption("-cembalo") then
			gridName = util.GetParam("-grid", "cembalo/cembalo_platte_2d_sci_1.ugx")
			requiredSubsets = {"Inner", "Boundary"}
			boundaries = {"Boundary"}
		elseif util.HasParamOption("-33k") then
			gridName = util.GetParam("-grid", "amgtests/unstructured_33k.ugx") -- -???
			requiredSubsets = {"Inner", "Left", "Right", "Top", "Bottom"}
			boundaries = {"Bottom"}
		elseif util.HasParamOption("-3k") then
			gridName = util.GetParam("-grid", "amgtests/unstructured_3k.ugx")
			requiredSubsets = {"Inner", "Left", "Right", "Top", "Bottom"}
			boundaries = {"Bottom"}
		elseif util.HasParamOption("-beam") then
			gridName = util.GetParam("-grid", "beamlong.ugx")
			requiredSubsets = {"Inner", "Boundary"}
			boundaries = {"Boundary"}
		else
			gridName = util.GetParam("-grid", "unit_square_01/unit_square_01_quads_2x2.ugx")
			requiredSubsets = {"Inner", "Boundary"}
			boundaries = {"Boundary"}
		end
		
		
		dom = MyLoadGrid({ dim=2, cpu = cpu, gridName = gridName, requiredSubsets = requiredSubsets} )
	end
	
		
	util.CheckAndPrintHelp("");
	
	
	if util.HasParamOption("-LE") then
		approxSpace, domainDisc = CreateLEDisc({Domain = dom, BoundaryName = boundaries[1]})		
	else
		approxSpace = CreateStdApproxSpace(dom)	
		domainDisc = CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries})
	end
else

	gridName = util.GetParam("-grid", "../skin/2D/geometry/brick/skin2d.aniso.ugx")
	requiredSubsets = {"LIP", "COR", "BOTTOM_SC", "TOP_SC"}
	dom = MyLoadGrid({ dim=2, gridName = gridName, requiredSubsets = requiredSubsets} )
	boundaries = {"TOP_SC", "BOTTOM_SC"}

	approxSpace = CreateStdApproxSpace(dom)
	domainDisc = CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries,
		diffusion={LIP=1, COR=10}
		})
end
OrderCuthillMcKee(approxSpace, true)
	

--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
print ("Setting up Algebra Solver")

convCheck = ConvCheck(5, 1e-5, 1e-5, false)
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
--SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")
print(" done.")

amg = FAMGPreconditioner()
SetFAMGCommandLineParameters(amg, {
	--delta=0.030473,
	delta=0.1,
	smooths=4,
	maxNodesForBase=2000	
	})


sgs = SymmetricGaussSeidel()
amg:set_presmoother(sgs)
amg:set_postsmoother(sgs)
--SetDebugLevel(debugID.AMG.Coarsening, 4)

if util.HasParamOption("-writeMatrices") then
	matrixWritePath=CreateOutputDirectory("RAMmatrices")
	--matrixWritePath = "/Users/mrupp/matrices/"
	amg:set_matrix_write_path(matrixWritePath)
	amg:set_write_f_values(true)
	--amg:set_debug_write_all_vectors(true)
end
amg:set_write_testvectors(true)

testVectorSmooths = util.GetParamNumber("-testVectorSmooths", 2)
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
amg:set_testvector_smooths(testVectorSmooths)
if util.HasParamOption("-tvsGS") then
	amg:set_testvector_smoother(SymmetricGaussSeidel()) --Jacobi(0.6))
else
	amg:set_testvector_smoother(Jacobi(0.6))
end

amg:set_delta(0.03)

amg:set_xc_calc_interpolation_routine(3)
--amg:set_delta(0.045)
amg:set_delta(0.2)
amg:set_delta(util.GetParamNumber("-delta", amg:get_delta()))
amg:set_deficit_solve_tolerance(1e-15)

amg:set_XC_LE_griebel(true)

amg:set_fsmoothing(true)
amg:set_delta(0.001, 4)
--amg:set_verbose(5)
--amg:add_watchnode(2, 939)
--amg:add_watchnode(1, 2800)
--amg:add_watchnode(0, 1084)
-- amg:add_watchnode(1, 839)
--amg = RSAMGPreconditioner()

--amg:add_watchnode(1, 21066)

--amg = RSAMGPreconditioner()
amg:set_strong_connection_external(0.1)
convCheck = ConvCheck(30, 1e-12, 1e-14)
linSolver = LinearSolver(amg, convCheck)
solver = linSolver

print(solver:config_string())
print("Init Solver...")

--SetDebugLevel(debugID.AMG, 5)
--SetDebugLevel(debugID.AMG.FAMG.Interpolation, 0)
amg:set_verbose(5)

if util.HasParamOption("-beam") then
	amg:set_beam(true)
end
solver:init(A)

lis=amg:level_information_str()
	print(lis)
print("Solve...")


--solve_random(solver, approxSpace)

amg:set_one_init(true)


if util.HasParamOption("-beam") then
	amg:set_fsmoothing(false)

	sm = SparseBlockGaussSeidel2(2)
	bgsRes = test_smoothsteps(sm, {1, 2, 4, 16})
	
	
else
	sgsRes = test_smoothsteps(SymmetricGaussSeidel(), {1, 2, 4, 8})
	--jacRes = test_smoothsteps(Jacobi(0.6), {1, 2, 4, 8})	
end

if true then
	res = solve_lin(LinearSolver(), Jacobi(0.6), "jac")
	print(res)
	stretchX = util.GetParamNumber("-stretchX", 1.0)
	stretchY = util.GetParamNumber("-stretchY", 1.0)
	print("scaling domain by ("..stretchX..", "..stretchY..")") 
	ScaleDomain(dom, MakeVec(0,0,0), MakeVec(1/stretchX, 1/stretchY, 1.0))
	SaveVectorForConnectionViewer(res.u, "smooth.vec")
	SaveVectorForConnectionViewer(res.d, "smooth_d.vec")
end


if ProcRank() == 0 then
	stats = {
	
	{ "procs", NumProcs() },
	--{ "numPreRefs", numPreRefs},
	{ "numRefs", numRefs },
	{ "dim", dim},
	{ "gridName", gridName},
	{ "main [ms]", GetProfileNode("main"):get_avg_total_time_ms()},
	{ "perform_refinement [ms]", GetProfileNode("perform_refinement"):get_avg_total_time_ms()},
	 
	{ "date", os.date("y%Ym%md%d") },
	{ "SVN Revision", GetSVNRevision()},
	{"host", GetBuildHostname()},
	{"commandline", util.GetCommandLine() },
	{"testVectorSmooths", testVectorSmooths} }
	
	table.append(stats, getAMGStats(amg)) 
	
	if sgsRes ~= nil then flattenTable(sgsRes, "sgs", stats) end
	if jacRes ~= nil then flattenTable(jacRes, "jac", stats) end
	if bgsRes ~= nil then flattenTable(bgsRes, "bgs", stats) end	
	
	flattenTable(getAMGLevelTable(amg), "level", stats)
	
	util.printStats(stats)
	if outfile ~= nil  then	
		util.writeFileStats(stats, outfile)
	end

end
	
WriteProfileData("pd.pdxml")

amg:check_vector(TestvectorsConstant(), "constant")
if dim == 2 then amg:check_vector(TestvectorsLinear2d(), "tvLin")
else amg:check_vector(TestvectorsLinear3d(), "tvLin") end


if true then
	amg:set_preiterations_mimum_defect_at_check(1e-9)
	amg:set_nr_of_preiterations_at_check(5)
	amg:set_num_postsmooth(4)
	amg:set_num_presmooth(4)
	--amg:set_postsmoother(SymmetricGaussSeidel())
	--amg:set_presmoother(SymmetricGaussSeidel())
	--amg:init_smoothers()
	if NumProcs() == 1 then amg:check() end
end
