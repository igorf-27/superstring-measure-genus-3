load("functions/theta-symbolic.sage")
load("functions/theta-numeric.sage")
load("functions/modular-forms.sage")
load("functions/modular-forms-auxiliary.sage")
load("functions/detrep-ex27.sage")
load("functions/invariants.sage")
load("functions/invariants-auxiliary.sage")
load("functions/invariants-from-files.sage")

prec = 200
CC = ComplexField(prec)

#tau = 1/5*matrix([[2+6*I,1-2*I,-1-3*I],[1-2*I,-2+4*I,-3+I],[-1-3*I,-3+I,3+4*I]]).change_ring(CC)
tau = random_period_matrix(3)
thetas = numeric_thetas(tau)

print(" * We pick at random a 3x3 period matrix $\\tau$ and compute $\\Xi_3(\\tau)$")
print(" * first by using Grushevsky's basis [Morozov08],")
print(" * then by using our formula (1.5b) in terms of invariant theory of nets")
print(" * with coeffiients $k_1 = 3/2, k_2 = 0, k_3 = -1/2$.")
print(Xi3(thetas,"grushevsky"))
print(Xi3(thetas,"invariants_of_nets"))
print("\n")

tau1 = random_period_matrix(1)
tau2 = random_period_matrix(2)
tau3 = block_diagonal_matrix(tau1,tau2)

thetas1 = numeric_thetas(tau1)
thetas2 = numeric_thetas(tau2)
thetas3 = numeric_thetas(tau3)

print(" * We pick at random a 1x1 period matrix $\\tau_1$ and a 2x2 one $\\tau_2$.")
print(" * We check that formula (1.1) for $\\Xi_1$ for from our paper")
print(" * and the formula in terms of Grushevsky's basis [Morozov08]")
print(" * describe the same modular form by evaluating both on $\\tau_1$.")
Xi1_sh = Xi1(thetas1,"short_formula")
print(Xi1_sh)
Xi1_gr = Xi1(thetas1,"grushevsky")
print(Xi1_gr)
print("\n")


print(" * Now we likewise compare formula (1.2) for $\\Xi_2$ from our paper")
print(" * with the formula in terms of Grushevsky's basis [Morozov08]")
print(" * by evaluating both on $\\tau_2$.")
Xi2_sh = Xi2(thetas2,"short_formula")
print(Xi2_sh)
Xi2_gr = Xi2(thetas2,"grushevsky")
print(Xi2_gr)
print("\n")


print(" * Now we let $\\tau_3$ be the 3x3 block-diagonal matrix with blocks $\\tau_1$")
print(" * and $\\tau_2$, and we check the factorization condition (1.3). That is,")
print(" * we check that $\\Xi_3(\\tau_3) = \\Xi_1(\\tau_1) \\Xi_2(\\tau_2)$.")
print(" * We calculate $\\Xi3$ via Grushevsky's basis.")
print(" * (Note that our formula (1.5b) is not suitable for block-diagonal")
print(" * matrices because it will try to calculate 0/0 in this case.)")
print(Xi3(thetas3,"grushevsky"))
print(Xi1_sh * Xi2_sh)

