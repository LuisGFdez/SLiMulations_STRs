
import numpy as np
import pandas as pd
from scipy.stats import geom
import pandas as pd
import sys
import argparse


def alpha_microsat(i, u, v, lo):
        return max(0, min(1, u - (v * (i - lo))))

def beta_microsat(i, s, mu ,lo):
        """
        beta helper function from Sainudiin et al. (2004)
        represents the mutation rate from allele i
        :param float s: strength of length dependence on mutation rate
        (i.e. s=0 uniform rate)
        :param int lo: lower bound
        :param int hi: upper bound
        """
        ##return mu*(1 + ((i - lo) * s))
        return 1 + ((i - lo) * s)
        
    
def gamma_microsat(m, i, j, lo, hi):
        """
        gamma helper function from Sainudiin et al. (2004)
        :param float m: prob success
        :param int i: in-state
        :param int j: out-state
        :param int lo: lower bound of microsat repeat length distribution
        :param int hi: upper bound of microsat repeat length distribution
        """
        num = m * (1 - m) ** (np.abs(i - j) - 1)
        if lo <= i < j <= hi:
            denom = 1 - (1 - m) ** (hi - i)
        elif lo <= j < i <= hi:
            denom = 1 - (1 - m) ** (i - lo)
        else:
            raise TypeError("microsat model error. i and j must be between lo and hi")

        return num / denom
 


hi = int(sys.argv[1]) 
lo = int(sys.argv[2]) 
folder=str(sys.argv[3])
##print(hi)
##print(lo)
num_alleles = hi - lo + 1
ancestral_alleles= np.arange(lo,hi+1,1)

##print (ancestral_alleles)
rate_matrix = np.zeros((num_alleles, num_alleles))
mutation=np.zeros((num_alleles, num_alleles))
##Parameter to simulate a stepwise mutation model
s=0
u=0.5
v=0
p=1
m=1
##u=0.68
##v=0.043
##m=0.01
##p=0.2
mu=1e-6

for i in range(lo, hi + 1):
    
    for j in range(lo, hi + 1):
            # shift coords to account for 0-indexing in rate matrix
        ii = i - lo
        jj = j - lo
        if i == j - 1:
            rate_matrix[ii, jj] = (
                  beta_microsat(i, s, mu ,lo)
                    * alpha_microsat(i, u, v, lo)
                    * (p + ((1 - p) * gamma_microsat(m, i, j, lo, hi)))
                )
            ##print ("mutation rate", beta_microsat(i, s, mu ,lo))
            #print("probability of expansion:", alpha_microsat(i, u, v, lo))
            #print("# Multistep units:", gamma_microsat(m, i, j, lo, hi) )
            #print(beta_microsat(i, s, mu ,lo) * alpha_microsat(i, u, v, lo) * (p + ((1 - p) * gamma_microsat(m, i, j, lo, hi))))
        elif i < j - 1:
                rate_matrix[ii, jj] = (
                    beta_microsat(i, s, mu ,lo)
                    * alpha_microsat(i, u, v, lo)
                    * (1 - p)
                    * gamma_microsat(m, i, j, lo, hi)
                )
                ##print ("mutation rate", beta_microsat(i, s, mu ,lo))
        elif i == j + 1:
               #print ("i:",i)
                #print ("j:",j)
                rate_matrix[ii, jj] = (
                    beta_microsat(i, s, mu ,lo)
                    * (1 - alpha_microsat(i, u, v, lo))
                    * (p + ((1 - p) * gamma_microsat(m, i, j, lo, hi)))
                )
                ##print ("mutation rate", beta_microsat(i, s, mu ,lo))
        elif i > j + 1:
                rate_matrix[ii, jj] = (
                    beta_microsat(i, s, mu ,lo)
                    * (1 - alpha_microsat(i, u, v, lo))
                    * (1 - p)
                    * gamma_microsat(m, i, j, lo, hi)
                )
                ##print ("mutation rate", beta_microsat(i, s, mu ,lo))
        else:
                rate_matrix[ii, jj] = 0
##

##with mutation rate
    # scale by max row sum; then normalize to 1
row_sums = rate_matrix.sum(axis=1, dtype="float64")

alpha = np.max(row_sums)
rate_matrix /= alpha
row_sums = rate_matrix.sum(axis=1, dtype="float64")
    # the max(0, *) is to avoid floating point error
np.fill_diagonal(rate_matrix, np.fmax(0.0, 1.0 - row_sums))

transition_matrix=rate_matrix
##transition_matrix=transition_matrix.flatten()
####root distribution
 # solve for stationary distribution

##S, U = np.linalg.eig(transition_matrix.T)
##U = np.real_if_close(U, tol=1)
##stationary = np.array(U[:, np.where(np.abs(S - 1.0) < 1e-8)[0][0]])
##stationary = stationary / np.sum(stationary)
##root_distribution = stationary

transition_matrix=pd.DataFrame(transition_matrix)
transition_matrix.to_csv(f'{folder}/transition_matrix.csv',header=False, index=False)
ancestral_alleles.tofile(f'{folder}/ancestral_alleles.csv',sep=",")
#root_distribution=pd.DataFrame(root_distribution)
##root_distribution.tofile('root_distribution.csv',sep=",")


##print(transition_matrix)
##print(root_distribution)

