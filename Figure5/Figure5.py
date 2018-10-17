import numpy as np
import sdeint
import pylab as plt
from scipy import integrate
from sklearn.preprocessing import normalize

####### This is the code to reproduce Figure 5 in the main text of the paper. 
####### The first half of the code compute the extinction probability over an ensemble of 
####### 250 realization (slow)
####### The secon half compute the time scale separation within the 2D feasibility domain

realization = 250
expected_transition= []; exploratory_distance = []
SS = 600

#############
d = 2
#############
### Few auxiliary functions
def distance(X,Y):
    return(np.sum(abs(X - Y)));
def separation(X):
   lambda_ = np.linalg.eigvals(X)
   return((max(lambda_) - min(lambda_))/max(lambda_));

#### Create the interaction matrix
#### Make sure that the feasibility domain is not too
#### large or too small and that the fixed point is stable
run = 0;
while run == 0:
    A = np.exp(np.random.rand(d,d))
    A = normalize(A, axis = 0, norm = 'l1')
    A_k = np.array(np.repeat(1, d))
    A_k = np.matrix(np.diag(A_k))
    A = A - np.diag(np.diag(A)) + A_k
    A = normalize(A, axis = 0, norm = 'l1')
    if np.linalg.det(A) > 0.4 and np.linalg.det(A) < 0.51:
        run = 1


A = np.asmatrix(A) ### Set it as matrix type for later in the SDE
starting_point = np.array([A[0,0], A[1,0]]) ### One column
arrival_point = np.array([A[0,1], A[1,1]])  ### The opposite column
max_dist = distance(starting_point, arrival_point)  ### Distance to travel
##################################
print '--- Dimensions?', d
print 'target distance:', max_dist
######################################################
######################################################
###### Generate panel (a) (slow)
n_steps = 80; #### Number of step to go from one side of the other of the domain
DeltaB = max_dist/n_steps; ### step size
extinction_probability = True
if extinction_probability == True:
    run = 0
    b_increment = 0;
    rlz = 0

    while run == 0:
        ext = np.zeros(realization)
        b = np.array([A[0,0] - b_increment, A[1,0] + b_increment])
        b = b/np.sum(b)
        t = np.linspace(0, 3000,  5000)
        for i in range(0,realization):    
            def f(X, t):
                dydt = np.array([X[s]*(b[s] - np.sum(np.dot(A,X)[0,s])) for s in range(0,len(X))])
                return dydt
            def u(X, t):
                dydt = np.array([1./np.sqrt(SS) * np.sqrt(X[s]*(b[s] + np.sum(np.dot(A,X)[0,s]))) for s in range(0,len(X))])
                return np.diag(dydt)
        
       
            x0 = np.array(np.repeat(1.,d))
            
            stochastic_result = sdeint.itoint(f, u, x0, t)
            for cl in range(0,d):
                clean_ =  np.nan_to_num(stochastic_result[:,cl])
                if min(clean_) <= 0:
                    ext[i] = 1
                    
        expected_transition.append(np.sum(ext)/len(ext))
        distance_from_border = distance(b, starting_point)
        exploratory_distance.append(distance_from_border)
        print distance_from_border, expected_transition[rlz]
        if distance_from_border >= max_dist:
            run = 1
        else:
            b_increment +=  DeltaB;
            rlz += 1;    
    nome_file = 'FD_exploration_in_%i' % d + 'd_large.txt' 
    s = open(nome_file, 'w')
    [s.write('%f %f\n' % (exploratory_distance[i], expected_transition[i])) for i in range(0,len(expected_transition))]
    s.close()
    fig = plt.figure()
    plt.plot(expected_transition, 'ro')
    plt.xlabel('Distance from left border')
    plt.ylabel('Extinction probability')

    
############################################
######## Second half of the code
######## Here you generate Figure 5 panel (b)
eg = []
exploratory_distance = []
b_increment = 0
run = 0;
while run == 0:
    b = np.array([A[0,0] - b_increment, A[1,0] + b_increment])
    b = b/np.sum(b)
    deterministic_result = np.array(np.dot(np.linalg.inv(A),b));
    distance_from_border = distance(b, starting_point)
    exploratory_distance.append(distance_from_border)
    abundances = np.diag(deterministic_result[0])
    J = np.matmul(abundances, A)
    eg.append(separation(J))
    if distance_from_border >= max_dist:
        run = 1
    else:
        b_increment +=  DeltaB;
nome_file = 'FD_exploration_time_scales_in_%i' % d + 'd_large.txt' 
s = open(nome_file, 'w')
[s.write('%f %f\n' % (exploratory_distance[i], eg[i])) for i in range(0,len(exploratory_distance))]
s.close()
fig2 = plt.figure()
plt.plot(exploratory_distance, eg, ls = '', marker = 'o')
plt.xlabel('Distance from left border')
plt.ylabel('Time Scale separation')

plt.show()
    


