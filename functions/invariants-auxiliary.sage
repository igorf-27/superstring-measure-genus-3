def terms_of_polynomial(f):
	c = f.coefficients()
	m = f.monomials()
	l = []
	for i in range(len(c)):
		l.append(c[i]*m[i])
	return l



def coefficients_of_quartic(q):
	R = q.parent()
	q_coefficients = {}
	for i in range(5):
		for j in range(5-i):
			k = 4-i-j
			mon = (R.0)^i * (R.1)^j * (R.2)^k
			coef = q.coefficient(mon).constant_coefficient()
			q_coefficients.update({(i,j,k):coef})
	return q_coefficients



def entries_of_symmetric_matrix(A):
	size = A.ncols()
	l = []
	for i in range(size):
		for j in range(i,size):
			l.append(A[i,j])
	return l



def entries_of_general_net(A,B,C):
	l = []
	for M in [A,B,C]:
		l.extend(entries_of_symmetric_matrix(M))
	return l



def entries_of_ex27_net(A,B,C):
	l = []
	for i in range(1,4):
		for M in [A,B,C]:
			l.append(M[0,i])
	return l



def quartic_from_net(A,B,C):
	R = A.base_ring()
	S.<x,y,z> = R[]
	q = (x*A+y*B+z*C).det()
	return q



def generic_symmetric_matrix(R, size, number_of_first_variable):
	A = matrix(R, size)
	n = number_of_first_variable
	for i in range(size):
		for j in range(size):
			if i<=j:
				A[i,j]=R.gens()[n]
				n=n+1
			else:
				A[i,j]=A[j,i]
	return A



# Takes an element poly of the Weyl algebra of some ring;
# poly must be a polynomial (i. e. a differential operator of order 0).
# Returns the element of the same Weyl algebra obtained from poly
# by changing $x$ into $\frac\partial{\partial x}$ for all variables.
def differential_operator(poly):
	W = poly.parent()
	R = W.polynomial_ring()
	n = len(R.gens())
	l = poly.list()
	d = 0
	for encoded_summand in l:
		summand = 1
		for i in range(n):
			summand = summand * (W.gen(n+i))^encoded_summand[0][0][i]
		summand = summand * encoded_summand[1]
		d = d + summand 
	return d

