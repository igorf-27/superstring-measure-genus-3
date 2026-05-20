# This generates files in ./auto-generated/invariants-as-polynomials/

load("functions/invariants.sage")
load("functions/invariants-auxiliary.sage")
load("functions/invariants-by-algorithms.sage")
load("functions/invariants-by-algorithms-auxiliary.sage")

# Loaded for 2 purposes:
# 1) to compute I3 for nets after I3 for quartics has been computed and written to a file, 
# 2) path_to_precomputed_invariants is defined there.
load("functions/invariants-from-files.sage")


print("Computing $I_3$ for quartics, this may be slow...")
l=[]
for i in range(5):
	for j in range(5-i):
		l.append('f'+str(i)+str(j)+str(4-i-j))
R = PolynomialRing(QQ,l)
S.<x,y,z> = R[]

q = 0
for i in range(5):
	for j in range(5-i):
		k = 4-i-j
		q = q + x^i*y^j*z^k*R('f'+str(i)+str(j)+str(k))

with open(path_to_precomputed_invariants + "I3-quartic.txt", "w") as f:
	f.write(str(I3_quartic(q,"differential"))+"\n")


print("Computing invariants for nets...")

R.<a,b,c,e,f,g,p,q,r> = QQ[]
A=matrix([[0,a,e,p],[a,0,0,0],[e,0,0,1],[p,0,1,0]])
B=matrix([[0,b,f,q],[b,0,0,1],[f,0,0,0],[q,1,0,0]])
C=matrix([[0,c,g,r],[c,0,1,0],[g,1,0,0],[r,0,0,0]])

with open(path_to_precomputed_invariants + "Lambda-ex27-net.txt", "w") as f:
	f.write(str(Lambda_net(A,B,C, "pfaffian"))+ "\n")
with open(path_to_precomputed_invariants + "I3-ex27-net.txt", "w") as f:
	f.write(str(I3_net(A,B,C, "quartic", algorithm_for_quartic="from_file"))+"\n")
with open(path_to_precomputed_invariants + "Qprime-ex27-net.txt", "w") as f:
	f.write(str(Qprime_net(A,B,C, "gizatullin"))+"\n")
with open(path_to_precomputed_invariants + "I-ex27-net.txt", "w") as f:
	f.write(str(I_net(A,B,C, "gizatullin").factor())+"\n")
with open(path_to_precomputed_invariants + "J-ex27-net.txt", "w") as f:
	f.write(str(J_net(A,B,C, "ex27_gizatullin").factor())+"\n")

