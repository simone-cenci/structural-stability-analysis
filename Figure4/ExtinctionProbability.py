import numpy as np
import sdeint
import pylab as plt
from scipy import integrate
from sklearn.preprocessing import normalize

############################################################
####### Run as python -W ignore ExtinctionProbability.py  ##
############################################################
############################################################
###### This file generate Figure 4 (main Figure)
###### Compute the transition probability as a function of the size of 
###### the feasibiltiy domain for a symmetric interaction matrix
###### Note that the code will be slow becuase of the integration of the
###### SDE and the large number of realizations needed in order to
###### get a clean curve for the extinction probability
#############################################################
expected_transition= []; FeasibilitySize = []
#### Off-diagonal interactions
alpha = 0.6
#### Number of realizations for each size of the feasibility domain
realization = 800
#### System size
SS = 300
#### Number of dimensions
d = 2
#### Choose if work with the normalized FD. In this case choose the norm
#### The norm is also needed to compute the size of the feasibility domain
Normalized = False
#### Set Norma and normalization together, i.e., 'lq', q
Norma, normalization = 'l1', 1
#### Functions to compute the size of the feasibility domain
def ColumnNormProduct(X, NormType):
    return(np.prod(np.linalg.norm(X,NormType,axis = 0)))
def SizeFeasibilityDomain(X, NormType):
    return(np.abs(np.linalg.det(X))/(ColumnNormProduct(X, NormType)))
###
k = 0.
data_points = 40
#### The time of integration as to be big enough to observe extinction. At least much greater than SS
t = np.linspace(0, 3000,  5000)
print 'Extinction Probability in', d, 'dimensions'
for z in range(0,data_points):
    ext = np.zeros(realization)
    b = np.array(np.repeat(1.,d))
    #### Create the interaction matrix
    A = np.matrix(np.full((d,d), alpha))
    A_k = np.array(np.repeat(k, d))
    A_k = np.matrix(np.diag(A_k))
    A = A + A_k
    ### If Normalize is True then the interaction matrix and the Growth rate are normalized
    if Normalized == True:
        b = b/np.linalg.norm(b, normalization)
        A =  normalize(A, axis = 0, norm = Norma)
        A = np.asmatrix(A)
    #######################################################
    def f(X, t):
        dydt = np.array([X[s]*(b[s] - np.sum(np.dot(A,X)[0,s])) for s in range(0,len(X))])
        return dydt
    def u(X, t):
        dydt = np.array([1./np.sqrt(SS) * np.sqrt(X[s]*(b[s] + np.sum(np.dot(A,X)[0,s]))) for s in range(0,len(X))])
        return np.diag(dydt)
        
       
    x0 = np.array(np.repeat(1.,d))
    for i in range(0,realization):
        stochastic_result = sdeint.itoint(f, u, x0, t)

        for cl in range(0,d):
            clean_ =  np.nan_to_num(stochastic_result[:,cl])
            if min(clean_) <= 0:
                ext[i] = 1
        
    expected_transition.append(np.sum(ext)/len(ext))
    FeasibilitySize.append(SizeFeasibilityDomain(A, normalization))
    k = k + 0.002
    print SizeFeasibilityDomain(A, normalization), np.sum(ext)/len(ext)
    
nome_file = 'TransitionProbability_in_%i' % d + 'd.txt' 
s = open(nome_file, 'w')
[s.write('%f %f\n' % (FeasibilitySize[i], expected_transition[i])) for i in range(0,len(expected_transition))]
s.close()
plt.plot(FeasibilitySize, expected_transition, 'ro')
plt.show()
    


