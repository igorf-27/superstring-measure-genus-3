# The Toeplitz invariant $\Lambda$ as a polynomial in 30 variables.
def Lambda_polynomial(R):
	S=PolynomialRing(R,'x',30)
	A=generic_symmetric_matrix(S,4,0)
	B=generic_symmetric_matrix(S,4,10)
	C=generic_symmetric_matrix(S,4,20)
	Z=zero_matrix(S,4)
	M=block_matrix([[Z,-C,B],[C,Z,-A],[-B,A,Z]],subdivide=False)
	return(M.pfaffian(algorithm="definition"))



# The contravariant net [Gizatullin07, Remark 2.6]
def contravariant_net(A,B,C):
	R=A.base_ring()
	t=Lambda_polynomial(R)
	S = t.parent()
	D=[matrix(R,4),matrix(R,4),matrix(R,4)]
	E=[matrix(S,4),matrix(S,4),matrix(S,4)]
	for n in range(3):
		G=generic_symmetric_matrix(S,4,10*n)
		for i in range(4):
			for j in range(4):
				E[n][i,j]=diff(t,G[i,j])
				if i!=j:
					E[n][i,j]=E[n][i,j] / 2
		for i in range(4):
			for j in range(4):
				D[n][i,j]=E[n][i,j](entries_of_general_net(A,B,C))	

	return D



# The four cubic forms [Gizatullin07, Section 4]
def cubics(A,B,C):
	R = A.base_ring()
	S.<u0,u1,u2,u3> = R[]
	quadrics = []
	for M in [A,B,C]:
		u_matrix = matrix([S.gens()])
		qf = (u_matrix * M * u_matrix.transpose())[0,0]
		quadrics.append(qf)
	J = matrix(S,3,4)
	for i in range(3):
		for l in range(4):
			J[i,l] = diff(quadrics[i],S.gen(l))
	cubics = []
	for k in range(4):
		columns = [r for r in range(4) if r != k]
		Jk = J[range(3),columns]
		cubics.append((-1)^k * Jk.det())
	return cubics



# Computes a certain 15x15 matrix whose entries are polynomials in 15 variables.
# The discriminant of a ternary quartic is the determinant of this matrix
# evaluated on the 15 coefficients of the quartic. See [Edge48].
def matrix_for_discriminant(q):
	R = q.base_ring()
	S.<a0,a1,a2,b0,b1,b2,c0,c1,c2,x0,x1,x2,y0,y1,y2>=R[]

	a = vector([a0,a1,a2])
	b = vector([b0,b1,b2])
	c = vector([c0,c1,c2])
	x = vector([x0,x1,x2])
	y = vector([y0,y1,y2])

	abc=matrix([a,b,c]).det()

	source_monomials = []
	for v in [a,b,c]:
		for i in range(5):
			for j in range(5-i):
				k = 4-i-j
				source_monomials.append(multinomial([i,j,k]) * v[0]^i * v[1]^j * v[2]^k)
	target_monomials = list(coefficients_of_quartic(q).values()) * 3

	T = 0
	for i in range(3):
		for j in range(3):
			k = 4-i-j
			if k<=2:
				T += (a*x)^i * (b*x)^j * (c*x)^k * (a*y)^(2-i) * (b*y)^(2-j) * (c*y)^(2-k) 
	T *=  (abc)^2 / 3
	terms_of_T = terms_of_polynomial(T)

	U = 0
	for i in range(len(terms_of_T)):
		u = terms_of_T[i]
		for j in range(len(source_monomials)):
			s = source_monomials[j]
			t = target_monomials[j]
			while u != 0 and s.divides(u):
				u = (u/s).numerator() * t
		U += u

	Rx.<x0,x1,x2> = R[]
	Rxy.<y0,y1,y2> = Rx[]
	U = U(0,0,0,0,0,0,0,0,0,x0,x1,x2,y0,y1,y2)

	entries_of_U = []	# we consider U as a symmetric bilinear form in y0, y1, y2
	for yy in [y0^2,y1^2,y2^2,y1*y2,y0*y2,y0*y1]:
		c = U.coefficient(yy).constant_coefficient() 	# to get from Rxy to Rx
		entries_of_U.append( c )
	for i in range(3,6):
		entries_of_U[i]  /=  2

	M = matrix(R,15)

	x = f.parent().gens()
	for i in range(3):
		for j in range(3):
			fij = x[j] * diff( f, x[i] )
			coefficients_of_fij = coefficients_of_quartic(fij)
			M[3*i+j]=list(coefficients_of_fij.values())	# dictionaries are ordered (since python v. 3.7)
	for i in range(6):
		w = entries_of_U[i]
		coefficients_of_w = coefficients_of_quartic(w) # we consider w as a quartic form in x0, x1, x2
		M[9+i]=list(coefficients_of_w.values())	# see above
	
	return M

