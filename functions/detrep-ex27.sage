# This file computes $A$ described in Appendix C of our paper.
# More precisely, $A$ is computed by the function detrep_ex27() below.

def theta_product(chars, thetas):
	return prod(theta_constant(ch,thetas) for ch in chars)

def adjugate_3x3(M):
	R = M.base_ring()
	N = matrix(R,3)
	for i in range(3):
		for j in range(3):
			A = M[[k for k in range(3) if k != i], [k for k in range(3) if k != j]]
			N[j,i] = A[0,0]*A[1,1] - A[0,1]*A[1,0]
			N[j,i] = N[j,i] * (-1)^(i+j)
	return N


# It will work for SR and CC only, otherwise most probably crash
def theta_ring(thetas):
	R_theta = thetas['th00'].base_ring()
	if R_theta != SR:
		R_theta = CC
	return R_theta
	

# The quantities $\beta_m$ of Appendix C.
def theta_gradiental_determinants(ch, thetas):
	T_ch = ['64','51','35']
	R_theta = theta_ring(thetas)
	T = matrix(R_theta,3)
	for i in range(3):
		for j in range(3):
			T[i,j] = theta_gradient(T_ch[i],j,thetas)
	T = adjugate_3x3(T)
	gr = vector([theta_gradient(ch,0,thetas), theta_gradient(ch,1,thetas), theta_gradient(ch,2,thetas)])
	return gr * T



def detrep_ex27(thetas):
	R_theta = theta_ring(thetas)
	S.<x0,x1,x2> = R_theta[]
	x_scaled = vector([ 
		theta_product(['04','40','67','76'],thetas) * theta_product(['03','12','24'],thetas) * x2,
		theta_product(['43','52','75'],thetas) * theta_product(['03','12','24'],thetas) * theta_product(['60'],thetas) * x1,
		theta_product(['43','52','75'],thetas) * theta_product(['04','40','67','76'],thetas) * x0])
	
	A=matrix(S,4)
	A[0,1] = theta_product(['04','41','50','66'],thetas) * theta_gradiental_determinants('77',thetas) * x_scaled
	A[0,2] = theta_product(['02','25','34','60'],thetas) * theta_gradiental_determinants('13',thetas) * x_scaled
	A[0,3] = theta_product(['01','04','10','37'],thetas) * theta_gradiental_determinants('26',thetas) * x_scaled
	A[1,2] = x2
	A[1,3] = x1
	A[2,3] = x0
	A += A.transpose()
	return [diff(A,x0).change_ring(R_theta),diff(A,x1).change_ring(R_theta),diff(A,x2).change_ring(R_theta)]


