import numpy as np
import os
from sklearn.preprocessing import normalize

#########################################################################
#### Run with python PlotSimplex.py 
#### It generate one panel of Figure 6 in the main text in a pdf file
#########################################################################
############ Functions needed
def separation(X):
   lambda_ = np.linalg.eigvals(X)
   return((max(np.real(lambda_)) - min(np.real(lambda_)))/max(np.real(lambda_)));
############################################################################
######## Create a stable interaction matrix and normalize it in the L-1 norm
d = 3
B = np.exp(np.random.rand(d,d))
B = normalize(B, axis = 0, norm = 'l1')
B_k = np.array(np.repeat(1, d))
B_k = np.matrix(np.diag(B_k))
B = B - np.diag(np.diag(B)) + B_k
B = normalize(B, axis = 0, norm = 'l1')
########
####### Just to get a sense of the interaction matrix
####### Note that the feasibility domain is just the determinant
####### Of the interaction matrix because the column are normalized in the L1 norm
print B
print np.linalg.eig(B)[0]
print '$\Omega$ = ', np.linalg.det(B), '--- Dimensions?', d
##################################################################
##### Initialize the vectors
b1 = []; b2 = []; b3 = []; eg = [];
##################################################################
for i in range(3000):
    #### Sample positive abundances and then take the corresponding growth rate
    x_ = np.random.uniform(1e-3,5,d)
    b = np.squeeze(np.asarray(np.dot(B,x_)))
    #### Normalize Growth rates
    b = b/np.sum(b)
    ##############################################
    ######### Compute the Jacobian
    abundances = np.matrix(np.diag(x_))
    J = np.matmul(abundances, B)
    ######### Compute the time scale separation
    eg.append(separation(J))
    b1.append(b[0]); b2.append(b[1]); b3.append(b[2]);

######### Write on file
s = open('Simplex.txt', 'w')
s.write('x y z density\n')
[s.write('%f %f %f %f\n' % (b1[i], b2[i], b3[i], eg[i])) for i in range(0,len(b1))]
s.close()
######### Create the Figure and clean the folder
os.system('pdflatex Figure.tex')
os.system('rm *.log *.aux')
os.system('evince Figure.pdf')
