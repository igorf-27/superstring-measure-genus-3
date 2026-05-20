path_to_precomputed_invariants = "auto-generated/invariants-as-polynomials/"



def read_file_as_string(file_name):
	with open(file_name) as file:
		str = file.read().replace('\n', ' ')
	return str



def I3_quartic_from_file(q):
	R = q.base_ring()
	for (i,j,k), coef in coefficients_of_quartic(q).items():
		exec("f" + str(i) + str(j) + str(k) + "=coef")
	s = read_file_as_string(path_to_precomputed_invariants + "I3-quartic.txt")
	return sage_eval(s,locals())



def invariant_ex27_net_from_file(A,B,C, invariant_name):
	a,b,c,e,f,g,p,q,r = entries_of_ex27_net(A,B,C)
	s = read_file_as_string(path_to_precomputed_invariants + invariant_name + "-ex27-net.txt")
	return sage_eval(s,locals())

