def I3_quartic_differential(q):
	R = q.parent().base_ring()
	S.<x0,x1,x2,y0,y1,y2,z0,z1,z2> = R[]
	q0 = q(x0,y0,z0)
	q1 = q(x1,y1,z1)
	q2 = q(x2,y2,z2)
	W = S.weyl_algebra()
	M = matrix([[x0,x1,x2],[y0,y1,y2],[z0,z1,z2]])
	d = differential_operator(W(M.det()^4))
	inv = R(d.diff(q0*q1*q2))
	return 1/(2^6 * 9) * inv



def discr_quartic_gordan(q):
	M=matrix_for_discriminant(q)
	return M.det()/2^21	



def Lambda_general_net_pfaffian(A,B,C):
	return Lambda_polynomial(A.base_ring())(entries_of_general_net(A,B,C))



def Qprime_general_net_gizatullin(A,B,C):
	R = A.base_ring()
	contranet = contravariant_net(A,B,C)
	cub = cubics(A,B,C)
	contracub_poly = cubics(contranet[0],contranet[1],contranet[2])
	W = cub[0].parent().weyl_algebra()
	contracub_diff = []
	for i in range(4):
		contracub_diff.append(differential_operator(W(contracub_poly[i])))

	sum_pkk = 0
	for k in range(4):
		sum_pkk = sum_pkk + R(contracub_diff[k].diff(cub[k]))
	return 1/2^7 * sum_pkk



def I_general_net_gizatullin(A,B,C):
	R=A.base_ring()
	S.<y0,y1,y2>=R[]
	T.<u0,u1,u2,u3>=S[]
	G=y0*A+y1*B+y2*C
	H=block_matrix( [ [G,matrix([[u0],[u1],[u2],[u3]])] ] ,subdivide=False)
	K=matrix([[u0,u1,u2,u3,0]])
	L=block_matrix( [ [H],[K] ] ,subdivide=False)
	d=L.det()
	M=matrix(R,10)
	u_monomials=[u0^2,u1^2,u2^2,u3^2,u0*u1,u0*u2,u0*u3,u1*u2,u1*u3,u2*u3]
	y_monomials=[y0^3,y1^3,y2^3,y0^2*y1,y0*y1^2,y1^2*y2,y1*y2^2,y2^2*y0,y2*y0^2,y0*y1*y2]
	for i in range(10):
		for j in range(10):
			M[i,j]=d.coefficient(u_monomials[i])(y0,y0,y0,y0).coefficient(y_monomials[j])(0,0,0)	#.change_ring(R)
	d=M.det()
	return 1/2^14 * d



def J_ex27_net_gizatullin(A,B,C):
	a,b,c,e,f,g,p,q,r = entries_of_ex27_net(A,B,C)
	d = matrix([[a,b,c],[e,f,g],[p,q,r]]).determinant()
	A = a*r^2-c*p*r
	E = a*f^2-b*e*f
	B = c*p*q + b*c*r - c*g*p - c*e*r + b*p*r + 2*a*g*r - 2*a*q*r - c^2*q
	G = b*c*f + c*e*f + b*e*g - b*f*p - b*e*q + 2*a*f*q - 2*a*f*g - b^2*g
	C = b*c*g - c*e*g + c*f*p + b*g*p + b*c*q + c*e*q - b*p*q + b*e*r - 2*a*g*q - 2*a*f*r + a*g^2 + a*q^2 - b^2*r - c^2*f
	F = b^2*r - c^2*f + b*c*g - b*c*q
	Jprime = matrix( [[4*A,3*B,2*C,G,0,0], [0,4*A,3*B,2*C,G,0], [0,0,4*A,3*B,2*C,G], [B,2*C,3*G,4*E,0,0], [0,B,2*C,3*G,4*E,0], [0,0,B,2*C,3*G,4*E]] ).determinant() / (16*F^2)
	J = (a*f*r*d)^2 * Jprime
	return J

