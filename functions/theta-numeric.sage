# This file is an interface to the package RiemannTheta.
from riemann_theta.riemann_theta import RiemannTheta

def numeric_thetas(tau):
	thetas = {}
	RTh=RiemannTheta(tau)
	g = tau.ncols()
	for ch in GF(2)^(2*g):
		lab = 'th' + char_human_readable(ch)
		if ch[:g]*ch[g:] == 0:
			thetas.update({ lab : RTh(char=ch) })
		else:
			thetas.update({ lab : 0 })
			for dir in range(g):
				lab2 = lab + str(dir)
				thetas.update({ lab2 : RTh(char=ch,derivs=[dir]) })
	return thetas

