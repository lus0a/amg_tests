
ug_load_script("../setup.lua")



InitSignals()

amgsetup.LoadStdDomAndApproxSpace2()


util.CheckAndPrintHelp("");


dom = MyLoadGrid({ gridName = gridName, requiredSubsets = requiredSubsets} )


if util.HasParamOption("-LE") then
	LinearElasticity = true
	approxSpace, domainDisc = CreateLEDisc({Domain = dom, BoundaryName = boundaries[1]})		
else
	approxSpace = CreateStdApproxSpace(dom)	
	domainDisc = CreateStdDisc({approxSpace=approxSpace, boundaries=boundaries})
end

A, u, b = AssembleAxb(approxSpace, domainDisc)
B = MatrixOperator()
v = GridFunction(approxSpace)
v:set(1.0)
domainDisc:adjust_solution(v)
domainDisc:assemble_mass_matrix(B, u)

amg = FAMGPreconditioner()


NSmooths = util.GetParamNumber("-nsmooths", 2)
SetFAMGCommandLineParameters(amg, {
	--delta=0.030473,
	delta=0.1,
	smooths=NSmooths,
	maxNodesForBase=2000,
	XC	
	})

amg:set_delta(util.GetParamNumber("-delta", 0.3))

amg:set_max_levels(2)
amg:set_max_nodes_for_base(1)
amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
amg:set_testvector_smooths(4)
amg:set_testvector_smoother(Jacobi(0.6))
amg:set_fsmoothing(false)
amg:set_deficit_solve_tolerance(1e-15)

if util.HasParamOption("-writeMatrices") then
	matrixWritePath=CreateOutputDirectory("RAMmatrices")
	amg:set_matrix_write_path(matrixWritePath)
	SaveMatrixForConnectionViewer(u, B, "mass.mat")
end
amg:set_write_testvectors(true)
	
smoother = SymmetricGaussSeidel()
amg:set_presmoother(smoother)
amg:set_postsmoother(smoother)
amg:set_num_presmooth(NSmooths)
amg:set_num_postsmooth(NSmooths)
	
	
--amg:set_testvectors(Test
if util.HasParamOption("-tvLin") then
	amg:set_testvectors(TestvectorsLinear2d())
end

if util.HasParamOption("-tvLE") then
	amg:set_testvectors(TestvectorsLinearElasticity2d())
end

-- solver
convCheck=ConvCheck(30, 1e-12, 1e-12, true) 
amgSolver = LinearSolver(amg, convCheck)
amgSolver:init(A)
amgSolver:apply(u, b)
print(amgSolver:config_string())

linSolver = amgSolver

-- kappa
function calcKappa2(tau, sigma)
	if tau <= 2*sigma then
		return tau/(4*sigma)
	else
		return 1-sigma/tau
		
	end
end

sigma = 1
tau = 0
bSmoother = true

---- \norm{Se}_1^2 \le \norm{e}_1^2 - \sigma \norm{e}_2^2

-- if off, this uses \norm{(I-PR')Se}_1^2 \le \tau \norm{e}_1^2
-- if on, without smoother
RSProof = false
--RSProof = true

globalSigma =1 
globalTau =0
-- this function approximates mu where
 -- \norm{(1-PR')e}_0^2 \le \tau \norm{e}_1^2
function calcTau(A, e, name)
	-- presmoothing
	Se = RandomFunction()
	Assign(Se, 1.0*e)
	 
	solver:apply_return_defect(Se, b:clone())
	
	D = MatrixDiagonal(A)
	DInv = MatrixDiagonalInverse(A)
	
	Se1 = VecProd(A*Se, Se)
	e0 = VecProd(D*e, e)
	e1 = VecProd(A*e, e)
	e2 = VecProd(DInv*(A*e), A*e)
	eNorm = VecProd(e, e)
	printAll=true
	if printAll then
		print("|Se|_1^2 = VecProd(A*Se, Se) = "..Se1)
		print("|e|^2 = VecProd(e, e) = "..eNorm)
		print("|e|_0^2 = VecProd(De, e) = "..e0)
		print("|e|_1^2 = VecProd(Ae, e) = "..e1)
		print("|e|_2^2 = VecProd(DInv*A*e, A*e) = "..e2)
	end
	
	if RSProof == true then
		IPRe = amg:check_interpolation(e, name..".vec")
		Se = e
	else
		IPRe = amg:check_interpolation(Se, name..".vec")
	end
	
	_e0 = VecProd(D*Se, Se)
	_e1 = VecProd(A*Se, Se)
	_e2 = VecProd(DInv*(A*Se), A*Se)
	_eNorm = VecProd(Se, Se)

	if printAll then
		print("|(I-PR')e|^2 = VecProd(e, e) = ".._eNorm)
		print("|(I-PR')e|_0^2 = VecProd(De, e) = ".._e0)
		print("|(I-PR')e|_1^2 = VecProd(Ae, e) = ".._e1)
		print("|(I-PR')e|_2^2 = VecProd(DInv*A*e, A*e) = ".._e2)
	end
	IPRe = _e0
	
	if printAll then
		print("IPRe = "..IPRe..", IPRe/e1 "..IPRe/e1)
	end
	
	local tau =0
	local sigma = 1
	if e1 > 1e-20 then
		tau = IPRe/e1
		--tau = math.sqrt(IPRe/e2)
	end 
	
	if Se1 > e0 and Se1 > 1e-12 then
		bSmoother = false
		print("---- NOT A SMOOTHER!  ---")
	end
	
	if e2 > 1e-20 then
		sigma = -(Se1-e1)/e2
	end
	
	globalTau = math.max(tau, globalTau)
	globalSigma = math.min(sigma, globalSigma)
	
	kappa2 = calcKappa2(tau, sigma)
	if printAll then
		print("---------------> tau here = "..tau .. " globalTau now "..globalTau)
		print("---------------> sigma here = "..sigma .. " globalSigma now "..globalSigma)
		print("---------------> kappa here = "..math.sqrt(kappa2))
	else
		write("tau = "..tau.." sigma = "..sigma.." ")
	end
end

function calcTau2(e, name)	
	calcTau(A, e, name)
end

convCheck1=ConvCheck(NSmooths, 1e-12, 1e-12, false) 
convCheck1:set_supress_unsuccessful(true)
solver=LinearSolver(SymmetricGaussSeidel(), convCheck1 )

solver:init(A)

-- test with random functions
print("Random test")
for i=1,10 do
	e = RandomFunction()
	domainDisc:adjust_solution(e)
	calcTau2(e, "random_"..i)
end

firstSigma=sigma
sigma =1

print("sin test")
N = u:size()
-- test with sin
if LinearElasticity then
	if N > 6 then N = 6 end
	for x1 = 0,N do
		for y1 = x1,x1 do
			e = CreateSinGridFunction2(x1, y1, x1, y1)
			calcTau(A, e, "__"..x1.."x1"..y1)
		end
	end
else
	if N > 20 then N = 20 end
	for x = 0,N do
		for y = x,x do --1,math.sqrt(u:size()) do
			e = CreateSinGridFunction(x, y)
			calcTau(A, e, "__"..x.."x"..y)
		end
	end
end


print("eigenvalue test")

numVecs = 15
eig = EigenSolver()

eig:set_linear_operator_A(A)
eig:set_linear_operator_B(B)
eig:set_max_iterations(20)
eig:set_precision(1e-12)
eig:set_preconditioner(ILUTScalar(0.0))
eig:set_pinvit(3)
eig:set_additional_eigenvectors_to_keep(4)
v={}
for k=1,numVecs do
	v[k] = RandomFunction()
	domainDisc:adjust_solution(v[k])
	eig:add_vector(v[k])
end

eig:apply()

for k=1,numVecs do
	SaveVectorForConnectionViewer(v[k], "ev_"..k..".vec")
	calcTau(A, v[k], "ev_"..k)
end

--end

if bSmoother == false then
	print("---- NOT A SMOOTHER!  ---")
end

function printKappa()
	print("\n\n")
	-- print kappa/tau
	kappa2 = calcKappa2(globalTau, globalSigma)
	kappa = math.sqrt(kappa2)
	print("\n\ntau = "..globalTau)
	print("\n\sigma = "..globalSigma)
	print("\n\nkappa2 = "..kappa2..", kappa = "..kappa)
	print("1-globalSigma/globalTau = "..1-globalSigma/globalTau)
	if tau <= 2*sigma then
		print("tau <= 2*sigma ")
		print("math.sqrt(tau/(4*sigma)) = "..math.sqrt(globalTau/(4*globalSigma)))
		print("math.sqrt(1-sigma/tau) = "..math.sqrt(math.abs(1-globalSigma/globalTau)))	
	else
	 	print("tau > 2*sigma ")
	 	print("math.sqrt(1-sigma/tau) = "..math.sqrt(math.abs(1-globalSigma/globalTau)))
	 	print("math.sqrt(tau/(4*sigma)) = "..math.sqrt(globalTau/(4*globalSigma)))
	end
	print("\n\n")
end

printKappa()

print(amg:config_string())
convCheck=ConvCheck(30, 1e-12, 1e-12, true)
amgSolver:set_convergence_check(convCheck)
resRandom=solve_random_d(amgSolver, approxSpace)
print(resRandom)

convCheck=ConvCheck(1, 1e-12, 1e-12, false)
convCheck:set_supress_unsuccessful(true)
amgSolver:set_convergence_check(convCheck)

-- get convergence in <Ae, e>-norm

local f = GridFunction(approxSpace)
f:set_random(-1, 1)		
local d = ZeroFunction()
nold = math.sqrt(VecProd(A*f, f))
AeConv = {}
for i=1,15 do
	amgSolver:apply(f, d)
	n = math.sqrt(VecProd(A*f, f))
	
	if n < 1e-12 then break end
	AeConv[i] = n/nold
	lastAeConv = n/nold
	
	nold =n
	ee = RandomFunction()
	Assign(ee, 1.0*f)
	calcTau(A, ee, "bla")
	--VecScaleAssign(f, 1/n, f)
	--VecScaleAssign(d, 1/n, d)
	--nold=1
end

for i=1,15 do
	print(AeConv[i])
end

printKappa()

print("\n")
print("kappa "..kappa)
print("lastConv Ae: "..lastAeConv)
print("lastReduction normal: "..resRandom.lastReduction)
print("\n\n\n")


print(amg:level_information_str())