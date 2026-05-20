def random_period_matrix(g):
	RR = RealField(prec)
	imag = matrix(CC,g)
	for i in range(g):
		imag[i,i]=RR.random_element(0,5)
	T = random_matrix(RR,g)
	imag = T * imag * T.transpose()
	real = random_matrix(RR,g) * 4
	real = real + real.transpose()
	return real + I*imag



def grushevsky_basis(g,thetas):
	if g>3:
		print("Error: Grushevsky basis for genus g>3 not implemented")
		return None

	zero_char = vector(GF(2),2*g)
	even_chars = [ch for ch in GF(2)^(2*g) if ch[:g]*ch[g:]==0]
	xi0 = theta_constant(zero_char,thetas)^16
	xi1 = theta_constant(zero_char,thetas)^8 * sum(theta_constant(ch,thetas)^8 for ch in even_chars)

	xi = [xi0,xi1]

	if g>=2:
		xi2 = theta_constant(zero_char,thetas)^4    *    sum(  theta_constant(ch1,thetas)^4   *   ( sum( theta_constant(ch2,thetas)^4 * theta_constant(ch1+ch2,thetas)^4  for  ch2 in even_chars) )   for ch1 in even_chars)
		xi.append(xi2)

	if g>=3:
		xi3 = theta_constant(zero_char,thetas)^2    *    sum(  theta_constant(ch1,thetas)^2   *    sum( theta_constant(ch2,thetas)^2 * sum( theta_constant(ch3,thetas)^2 * theta_constant(ch1+ch2,thetas)^2 * theta_constant(ch1+ch3,thetas)^2 * theta_constant(ch2+ch3,thetas)^2 * theta_constant(ch1+ch2+ch3,thetas)^2 for ch3 in even_chars) for ch2 in even_chars) for ch1 in even_chars)
		xi.append(xi3)

	return xi

