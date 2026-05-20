load("functions/invariants.sage")
load("functions/invariants-auxiliary.sage")
load("functions/invariants-by-algorithms.sage")
load("functions/invariants-by-algorithms-auxiliary.sage")
load("functions/invariants-from-files.sage")



def ex27(a,b,c,e,f,g,p,q,r):
	A = matrix([[0,a,e,p],[a,0,0,0],[e,0,0,1],[p,0,1,0]])
	B = matrix([[0,b,f,q],[b,0,0,1],[f,0,0,0],[q,1,0,0]])
	C = matrix([[0,c,g,r],[c,0,1,0],[g,1,0,0],[r,0,0,0]])
	return [A,B,C]



R.<x,y,z>=QQ[]
f = x^4+y^4+z^4
print(" * For the quartic ", f)
print(" * the invariant $I_3$ should be 144 and the discriminant 1.")
print("I_3 = ",I3_quartic(f,"differential"))	#144
print("discr = ", discr_quartic(f,"gordan"))	#1
print("\n")

f = x^4 - y^4 + 2*z^4 - 5*x^2*y*z
print(" * For the quartic ", f)
print(" * the invariant $I_3$ should be -288 and the discriminant -321499206081/524288.")
print("I_3 = ", I3_quartic(f,"differential"))	#-288
print("discr = ",discr_quartic(f,"gordan"))	#-321499206081/524288
print("\n")



R.<a,b,c,e,f,g,p,q,r>=QQ[]
parameters_ex51 = [r,b,c,b,r,g,c,g,r]
parameters_riemann = [-1-e-p,b,c,e,-1-b-q,g,p,q,-1-c-g]


print(" * Let us compute $\\Lambda, I, I_3$ and $Q'$ for the matrices")
print(" * of Example 5.1 in [Gizatullin07].\n")
print(" * $\\Lambda$ must be 0, see Example 5.1 (page 51).")
print("\\Lambda = ", Lambda_net(*ex27(*parameters_ex51),"ex27_from_file"))
print("\\Lambda = ", Lambda_net(*ex27(*parameters_ex51),"pfaffian"), "\n")

print(" * I is given by the formula of [Gizatullin07, Remark 7.7] (page 56).")
print("I = ", I_net(*ex27(*parameters_ex51),"ex27_from_file").factor())
print("I = ", I_net(*ex27(*parameters_ex51),"gizatullin").factor(),"\n")

print("$I_3$ must be proportional to the last formula of [Gizatullin07, Example 5.1] (p. 51):")
print("I_3 = ",I3_net(*ex27(*parameters_ex51),"ex27_from_file",algorithm_for_quartic="from_file").factor(),"\n")

print("$Q'$ must be 0, see [Gizatullin07, Example 5.1] (p. 51).")
print("Q' =", Qprime_net(*ex27(*parameters_ex51),"ex27_from_file"))
print("Q' =", Qprime_net(*ex27(*parameters_ex51),"gizatullin"), "\n")


print(" * Let us compute $J$ for the *Rieman normalization* of Example 2.7 (page 43).")
print(" * There must be 6 irreducible polynomials in its factorization, the first")
print(" * 5 polynomials must all have exponent 2, and the last one must have")
print(" * exponent 1. This is described in [Gizatullin07, Remark 10.2] (p. 61).")
print(" * The last polynomial is large, so we don't try to print it. We print the")
print(" * number of irreducible polynomials (must be 6), and then we try to print the")
print(" * 5 small polynomials with their exponents (technically, we ask Sage to print")
print(" * only those irreducible factors that consist of no more than 16 non-zero")
print(" * monomials).")

J_riemann = J_net(*ex27(*parameters_riemann),"ex27_gizatullin").factor()
n = len(J_riemann)
print("Number of irreducible polynomials: ",n)
s='J = '
for i in range(n):
	if len((J_riemann[i][0]).coefficients()) <= 16:
		s = s + '(' + str(J_riemann[i][0]) + ')' + '^' + str(J_riemann[i][1]) + ' * '
s = s + '...'
print(s)

