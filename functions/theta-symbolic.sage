def theta_constant(ch, thetas):
	if type(ch) != type(''):
		ch = char_human_readable(ch)
	return thetas['th'+ch]

def theta_gradient(ch, dir, thetas):
	if type(ch) != type(''):
		ch = char_human_readable(ch)
	return thetas['th'+ch+str(dir)]

def char_human_readable(ch):
	g = ch.length() / 2
	label_demi_length = len(str(2^g-1))
	a = 0
	b = 0
	for i in range(g):
		a += ZZ(ch[i]) * 2^(g-1-i)
		b += ZZ(ch[g+i]) * 2^(g-1-i)
	label_a = str(a)
	while len(label_a) < label_demi_length:
		label_a = '0' + label_a
	label_b = str(b)
	while len(label_b) < label_demi_length:
		label_b = '0' + label_b
	return label_a + label_b


# Takes a string of digits encoding a theta characteristic.
# Returns the corresponding 6-component vector with 0 or 1 as components.
# E. g. '13' gets convected into (0,0,1,0,1,1).
def char_binary(lab,g):
	label_demi_length = len(str(2^g-1))
	a = ZZ(lab[:label_demi_length])
	b = ZZ(lab[label_demi_length:])
	ch = vector(ZZ, 2*g)
	for i in range(g):
		a,r = a.quo_rem(2)
		ch[g-1-i] = r
		b,r = b.quo_rem(2)
		ch[2*g-1-i] = r
	return ch.change_ring(GF(2))


def symbolic_thetas(g):
	thetas = {}
	for ch in GF(2)^(2*g):
		lab = 'th' + char_human_readable(ch)
		if ch[:g]*ch[g:] == 0:
			thetas.update({lab:SR(lab)})
		else:
			thetas.update({lab:0})
			for dir in range(g):
				lab2 = lab + str(dir)
				thetas.update({lab2:SR(lab2)})
	return thetas

