
ug_load_script("laplace_sin_filtering.lua")

NSmooths = util.GetParamNumber("-nsmooths", 1)

if true then
	amg = FAMGPreconditioner()
	
	amg = createAMG1()
	amg:set_max_levels(2)
	amg:set_max_nodes_for_base(1)
	--amg:set_testvectors(TestvectorsExternal({ CreateSinGridFunction(1, 1) } ))
	
	
	smoother = Jacobi(0.6)
	smoother = SymmetricGaussSeidel()
	amg:set_presmoother(smoother)
	amg:set_postsmoother(smoother)
	if util.HasParamOption("-RSproof") then
		amg:set_num_presmooth(0)
	else
		amg:set_num_presmooth(NSmooths)
	end
	amg:set_num_postsmooth(NSmooths)
	amg:set_fsmoothing(true)
	
	convCheck=ConvCheck(30, 1e-12, 1e-12, true) 
	amgSolver = LinearSolver(amg, convCheck)
	amgSolver:init(A)
	amgSolver:apply(u, b)
	print(amgSolver:config_string())
	
	createLS(amg)
	
else
	amg = FAMGPreconditioner()
	
	
	NSmooths = util.GetParamNumber("-nsmooths", 2)
	SetFAMGCommandLineParameters(amg, {
		--delta=0.030473,
		delta=0.1,
		smooths=NSmooths,
		maxNodesForBase=2000,
		XC	
		})
	
	amg:set_max_levels(2)
	amg:set_max_nodes_for_base(1)
	amg:set_damping_for_smoother_in_interpolation_calculation(0.6)
	amg:set_testvector_smooths(4)
	amg:set_testvector_smoother(Jacobi(0.6))
	amg:set_fsmoothing(true)
	amg:set_deficit_solve_tolerance(1e-15)
	
	if util.HasParamOption("-writeMatrices") then
			matrixWritePath=CreateOutputDirectory("RAMmatrices")
			amg:set_matrix_write_path(matrixWritePath)
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
	
	
	-- solver
	convCheck=ConvCheck(30, 1e-12, 1e-12, true) 
	amgSolver = LinearSolver(amg, convCheck)
	amgSolver:init(A)
	amgSolver:apply(u, b)
	print(amgSolver:config_string())
	
	linSolver = amgSolver

end

function calcKappa2(tau, sigma)
	if RSProof or tau > 2*sigma then
		return 1-sigma/tau
	else
		return tau/(4*sigma)		
		
	end
end

sigma = 1
tau = 0
bSmoother = true

--- sigma: smoother constant. get minimum over all configurations (gets minimal for random vectors)
--- tau: interpolation constant. get maximum over all configurations. get max for smooth vectors

---- \norm{Se}_1^2 \le \norm{e}_1^2 - \sigma \norm{e}_2^2

-- if off, this uses \norm{(I-PR')Se}_1^2 \le \tau \norm{e}_1^2
RSProof = util.HasParamOption("-RSproof")
printAll = util.HasParamOption("-printAll")

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
	IPRe = _e0

	if printAll then
		print("|(I-PR')e|^2 = VecProd(e, e) = ".._eNorm)
		print("|(I-PR')e|_0^2 = VecProd(De, e) = ".._e0)
		print("|(I-PR')e|_1^2 = VecProd(Ae, e) = ".._e1)
		print("|(I-PR')e|_2^2 = VecProd(DInv*A*e, A*e) = ".._e2)
	end
	
	
	if printAll then
		print("IPRe = "..IPRe..", IPRe/e1 "..IPRe/e1)
	end
	
	local tau =0
	local sigma = 1
	---- interpolation property
	if e1 > 1e-20 then
		-- |(1-PR)e |^2_0 <= tau |e|^2_1   (RSAMG)
		-- |(1-PR)Se |^2_0 <= tau |e|^2_1  (FAMG)
		tau = IPRe/e1
		--tau = math.sqrt(IPRe/e2)
	end 
	
	if Se1 > e0 and Se1 > 1e-12 then
		bSmoother = false
		print("---- NOT A SMOOTHER!  ---")
	end
	
	---- smoothing property
	if e2 > 1e-20 then
		-- |Se|^2_1 <= |e|^2_1 - sigma |e|^2_2
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


function solve_sin2(x, y)
	local f = CreateSinGridFunction(x, y)
	domainDisc:adjust_solution(f)
	solver:apply_return_defect(f, b:clone())
	local avgRate = convCheck:avg_rate()
	if avgRate < 1e-12 then
		convMsg = avgRate.." - EXACT!"
	else
		convMsg = avgRate
	end	
	lastReduction = convCheck:defect()/convCheck:previous_defect()
	print("solve nuX = "..x..", nuY = "..nuY.." avgConv = "..convMsg.." lastReduction = ".. lastReduction )
	if false and avgRate < 1e-12 then return "EXACT", "EXACT", f
	else 
	return string.format("%2.2e", avgRate), string.format("%2.2e", lastReduction), f 
	end
end

convCheck=ConvCheck(2, 1e-12, 1e-12, false)
convCheck:set_supress_unsuccessful(true)


--solver=LinearSolver(GaussSeidel(), convCheck )
--solver=LinearSolver(BlockGaussSeidel(1), convCheck )
--solver=LinearSolver(ILUT(1e-7), convCheck )
--solver=LinearSolver(Jacobi(0.66), convCheck )
convCheck1=ConvCheck(NSmooths, 1e-12, 1e-12, false) 
convCheck1:set_supress_unsuccessful(true)
--smoother = util.GetParamFromList("-amgSmoother", "sgs", amgsetup.getSmootherList() or p.smoother)
solver=LinearSolver(smoother, convCheck1 )
solver:init(A)

print("\n\nRandom test\n---------------------")
e = RandomFunction()
DInv = MatrixDiagonalInverse(A)
for i=1,20 do
	domainDisc:adjust_solution(e)
	e = DInv*(A*e)
	Assign(e, (1/VecNorm(e)) *e)
end
domainDisc:adjust_solution(e)

for i=1,10 do
	e = RandomFunction()
	domainDisc:adjust_solution(e)
	calcTau2(e, "v")
end

firstSigma=sigma
sigma =1

print("\n\nsin test\n---------------------")
N = math.sqrt(u:size())
if N > 20 then N = 20 end
for x = 0,N do
	for y = x,x do --1,math.sqrt(u:size()) do
		e = CreateSinGridFunction(x, y)
		calcTau(A, e, "__"..x.."x"..y)
	end
end
--end

if bSmoother == false then
	print("---- NOT A SMOOTHER!  ---")
end

print("\n\nConfig\n---------------------")
print(amg:config_string())

print("\n\nLinear Solve\n---------------------")
convCheck=ConvCheck(30, 1e-12, 1e-12, true)
amgSolver:set_convergence_check(convCheck)
resRandom=solve_random_d(amgSolver, approxSpace, convCheck)
print(resRandom)

convCheck=ConvCheck(1, 1e-12, 1e-12, false)
convCheck:set_supress_unsuccessful(true)
amgSolver:set_convergence_check(convCheck)


local f = GridFunction(approxSpace)
f:set_random(-1, 1)		
local d = ZeroFunction()
nold = math.sqrt(VecProd(A*f, f))

if true then
	print("\n\nchecking remaining error\n---------------------")
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
	end
	
	for i, v in ipairs(AeConv) do
		print(v)
	end
end

print("\n\nCalculated Values\n---------------------")
if RSProof then
	print(" ---- RS PROOF ---- ")
else
	print(" ---- FAMG PROOF ---- ")
end
kappa2 = calcKappa2(globalTau, globalSigma)
kappa = math.sqrt(kappa2)
print("\n\ntau = "..globalTau)
print("\n\sigma = "..globalSigma)
print("\n\nkappa2 = "..kappa2..", kappa = "..kappa.."\n")
if tau <= 2*sigma then
	print("tau <= 2*sigma ")
	print("math.sqrt(tau/(4*sigma)) = "..math.sqrt(globalTau/(4*globalSigma)))
	print("math.sqrt(1-sigma/tau) = "..math.sqrt(math.abs(1-globalSigma/globalTau)))	
else
 	print("tau > 2*sigma ")
 	print("math.sqrt(1-sigma/tau) = "..math.sqrt(math.abs(1-globalSigma/globalTau)))
 	print("math.sqrt(tau/(4*sigma)) = "..math.sqrt(globalTau/(4*globalSigma)))
end
print("\n")
print("kappa "..kappa)
if lastAeConv then print("lastConv Ae: "..lastAeConv) end
print("lastReduction normal: "..resRandom.lastReduction)
print("\n\n\n")

print(amg:level_information_str())