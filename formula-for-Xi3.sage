# This prints our formula for $Xi_3$ in terms of theta constants and theta gradients
# obtained from invariant theory of nets. See Remark 4 in Section 1.7 of our paper.

load("functions/invariants.sage")
load("functions/invariants-auxiliary.sage")
load("functions/invariants-from-files.sage")
load("functions/theta-symbolic.sage")
load("functions/detrep-ex27.sage")
load("functions/modular-forms.sage")


thetas = symbolic_thetas(3)
print(Xi3(thetas,"invariants_of_nets"))

