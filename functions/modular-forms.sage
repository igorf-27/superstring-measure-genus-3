def chi18(thetas):
	return prod([theta_constant(ch,thetas) for ch in GF(2)^6 if ch[:3]*ch[3:] == 0])



def Xi1(thetas, algorithm="short_formula"):

	# Formula (1.1) of our paper.
	if algorithm == "short_formula":
		return theta_constant('00',thetas)^8 * theta_constant('01',thetas)^4 * theta_constant('10',thetas)^4

	# See [Morozov08].
	if algorithm == "grushevsky":
		xi0, xi1 = grushevsky_basis(1,thetas)
		return xi0 - 1/2*xi1



def Xi2(thetas, algorithm="short_formula"):

	# Formula (1.2) of our paper.
	if algorithm == "short_formula":
		return theta_constant('00',thetas)^4 * ( (theta_constant('03',thetas)*theta_constant('10',thetas)*theta_constant('21',thetas))^4 + (theta_constant('01',thetas)*theta_constant('12',thetas)*theta_constant('30',thetas))^4 + (theta_constant('02',thetas)*theta_constant('20',thetas)*theta_constant('33',thetas))^4 )
	
	# See [Morozov08].
	if algorithm == "grushevsky":
		xi0, xi1, xi2 = grushevsky_basis(2,thetas)
		return 2/3*xi0 - 1/2*xi1 + 1/12*xi2



def Xi3(thetas, algorithm="grushevsky"):

	# See [Morozov08].
	if algorithm == "grushevsky":
		xi0, xi1, xi2, xi3 = grushevsky_basis(3,thetas)
		return 8/21*xi0 - 1/3*xi1 + 1/12*xi2 - 1/168*xi3

	# Formula (1.5b) of our paper with k_1 = 3/2, k_2 = 0, k_3 = -1/2.
	if algorithm == "invariants_of_nets":
		A,B,C = detrep_ex27(thetas)
		Lambda = Lambda_net(A,B,C)
		I = I_net(A,B,C)
		J = J_net(A,B,C)
		Qprime = Qprime_net(A,B,C)
		return (3/2 * Lambda^3 - 1/2 * Qprime) * I / J * theta_constant('00',thetas)^16

