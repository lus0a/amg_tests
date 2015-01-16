
ug_load_script("laplace_sin_filtering.lua")

--amg = FAMGPreconditioner()
amg:set_injection_restriction(true)
amg:set_max_levels(2)
amg:set_max_nodes_for_base(1)
amg:set_testvectors(TestvectorsExternal({ CreateSinGridFunction(4, 4) } )) 
amgSolver = LinearSolver(amg, convCheck)
amgSolver:init(A)
amgSolver:apply(u, b)
print(amgSolver:config_string())
sigma = 1
tau = 0
bSmoother = true
function calcMu(A, Se, e)
	D = MatrixDiagonal(A)
	DInv = MatrixDiagonalInverse(A)
	
	Se1 = VecProd(A*Se, Se)
	e0 = VecProd(D*e, e)
	e1 = VecProd(A*e, e)
	e2 = VecProd(DInv*(A*e), A*e)
	
	print("|Se|_1^2 = VecProd(A*Se, Se) = "..Se1)
	print("|e|_0^2 = VecProd(De, e) = "..e0)
	print("|e|_1^2 = VecProd(Ae, e) = "..e1)
	print("|e|_2^2 = VecProd(DInv*A*e, A*e) = "..e2)
	
	IPRe = amg:check_interpolation(e, "v.vec")
	
	if e1 > 1e-12 and IPRe > tau*e1 then
		tau = IPRe/e1
	end 
	print("tau here: "..IPRe/e1)
	
	if Se1 > e0 and Se1 > 1e-12 then
		bSmoother = false
		print("---- NOT A SMOOTHER!  ---")
	end
	
	if Se1 > e0 - sigma*e2 and e0-Se1 > 1e-12 then
		sigma = (e0-Se1)/e2		 
	end
	
	print("---------------> tau = "..tau)
	print("---------------> sigma = "..sigma)
end

function calcMu2(e)
	Se = RandomFunction()
	Assign(Se, 1.0*e)
	 
	solver:apply_return_defect(Se, b:clone())
	calcMu(A, Se, e)
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

convCheck=ConvCheck(1, 1e-12, 1e-12, false)
convCheck:set_supress_unsuccessful(true)


solver=LinearSolver(GaussSeidel(), convCheck )
solver=LinearSolver(BlockGaussSeidel(1), convCheck )
solver=LinearSolver(ILUT(1e-7), convCheck )
solver=LinearSolver(Jacobi(0.66), convCheck )
--solver=LinearSolver(SymmetricGaussSeidel(), convCheck )

solver:init(A)

e = RandomFunction()
DInv = MatrixDiagonalInverse(A)
for i=1,20 do
domainDisc:adjust_solution(e)
e = DInv*(A*e)
Assign(e, (1/VecNorm(e)) *e)
end
domainDisc:adjust_solution(e)
calcMu2(e)

firstSigma=sigma
sigma =1

for x = 1,math.sqrt(u:size()) do
	for y = x,x do --1,math.sqrt(u:size()) do
		--y = x
		avg,last,Se = solve_sin2(x,y)
		e = CreateSinGridFunction(x, y)
		SaveVectorForConnectionViewer(e, x.."x"..y..".vec")
		
		
		calcMu(A, Se, e)
		
		SaveVectorForConnectionViewer(Se, "f_"..x.."x"..y..".vec")
		
		print(x .. ","..y.." = " ..last)
	end
end
--end

if bSmoother == false then
	print("---- NOT A SMOOTHER!  ---")
end
print("first sigma="..firstSigma.." 2nd sigma="..sigma)
print("---------------> tau = "..tau)

if tau <= 2*sigma then
	kappa2 = tau/(4*sigma)
else
	kappa2 = 1-sigma/tau
end
kappa = math.sqrt(kappa2)
print("kappa2 = "..kappa2..", kappa = "..kappa)

