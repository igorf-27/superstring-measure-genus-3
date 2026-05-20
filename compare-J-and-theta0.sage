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
print(" * We pick a 3x3 period matrix $\\tau$ at random and check whether")
print(" * relation (4.29) from our paper holds.")
tau = random_period_matrix(3)
thetas = numeric_thetas(tau)

A,B,C = detrep_ex27(thetas)
I = I_net(A,B,C)
J = J_net(A,B,C)
print("(I^8 / J^5)(A(\\tau))                   =", I^8 / J^5)
chi = chi18(thetas)
th0 = theta_constant('00',thetas)
print("(\\chi_{18}^2 / \\theta_{00}^{72})(\\tau) =", chi^2 / th0^72)

