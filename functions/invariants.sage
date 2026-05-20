# INVARIANTS OF TERNARY QUARTICS

def I3_quartic(q, algorithm = "from_file"):
	if algorithm == "differential":
		return I3_quartic_differential(q)
	if algorithm == "from_file":
		return I3_quartic_from_file(q)


# Macaulay2's discriminant is 2^(19+21) times bigger.
def discr_quartic(q, algorithm = "gordan"):
	if algorithm == "gordan":
		return discr_quartic_gordan(q)



def quartic_from_net(A,B,C):
	R = A.base_ring()
	S.<x,y,z> = R[]
	q = (x*A+y*B+z*C).det()
	return q



# INVARIANTS OF NETS
	
# The invariant $\Lambda$ [Gizatullin07, section 2]
def Lambda_net(A,B,C, algorithm = "ex27_from_file"):
	if algorithm == "pfaffian":
		return Lambda_general_net_pfaffian(A,B,C)
	if algorithm == "ex27_from_file":
		return invariant_ex27_net_from_file(A,B,C, "Lambda")


def I3_net(A,B,C, algorithm = "ex27_from_file", algorithm_for_quartic = "from_file"):
	if algorithm == "quartic":
		q = quartic_from_net(A,B,C)
		return 1/4 * I3_quartic(q, algorithm_for_quartic)
	if algorithm == "ex27_from_file":
		return invariant_ex27_net_from_file(A,B,C, "I3")


def Qprime_net(A,B,C, algorithm = "ex27_from_file"):
	if algorithm == "gizatullin":
		return Qprime_general_net_gizatullin(A,B,C)
	if algorithm == "ex27_from_file":
		return invariant_ex27_net_from_file(A,B,C, "Qprime")


# The Salmon invariant I [Gizatullin07, section 7]
def I_net(A,B,C, algorithm = "ex27_from_file"):
	if algorithm == "gizatullin":
		return I_general_net_gizatullin(A,B,C)
	if algorithm == "ex27_from_file":
		return invariant_ex27_net_from_file(A,B,C, "I")


# The tact invariant J [Gizatullin07, section 10]
def J_net(A,B,C, algorithm = "ex27_from_file", algorithm_for_discr = "quartic", algorithm_for_I = "gizatullin"):
	if algorithm == "discriminant":
		discr = discr_net(A,B,C, algorithm_for_discriminant)
		I = I_net(A,B,C, algorithm_for_I)
		return discr / I^2
	if algorithm == "ex27_gizatullin":
		return J_ex27_net_gizatullin(A,B,C)
	if algorithm == "ex27_from_file":
		return invariant_ex27_net_from_file(A,B,C, "J")

def discr_net(A,B,C, algorithm = "ex27_IIJ", algorithm_for_quartic = "gordan", algorithm_for_I = "ex27_from_file", algorithm_for_J = "ex27_from_file"):
	if algorithm == "quartic":
		q = quartic_from_net(A,B,C)
		return discr_quartic(q, algorithm_for_quartic)
	if algorithm == "ex27_IIJ":
		I = I_net(A,B,C, algorithm_for_I)
		J = J_net(A,B,C, algorithm_for_J)
		return I*I*J

