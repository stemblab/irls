#!puzlet

import numpy as np
import numpy.linalg as la
import time

# $\epsilon$-regularized IRLS solution to $\Phi u =b$ 
# http://stemblab.github.io/irls
 
# $\LaTeX$ definitions
# brackets: $$\def\lr#1{\left(#1\right)}$$

def solve(Phi,b,p=0,thresh=0.1,max_eps=1e-8,max_iter=200,eps_scale=10,
    eps_thresh=100):
   
    # initialization   

    epsilon=1 # $\epsilon$
    ul=la.lstsq(Phi,b)[0] # "last $u$", $u^{(n-1)}$ = LMS soln $\Phi u=b$
    m=ul.shape[0]
    Qn=np.zeros((m,m)) # $Q_n$

    # iteration

    for iter in np.arange(1,max_iter):

        # stop if error is small
        if epsilon<max_eps: break
    
        # main algorithm

        # $w_i=\lr{\lr{u_i^{(n-1)}}^2+\epsilon}^{p/2-1}$ 
    
    
        # $Q_i=1/w_i$
    
        np.fill_diagonal(Qn, 1./((np.abs(ul)**2+epsilon)**(p/2-1)))

        # $u^{(n)}=Q_n\Phi^H \lr{\Phi Q_n \Phi^H}^{-1}b$ 


        QPH = np.dot(Qn, Phi.conj().T) # $Q_n \Phi^H$
        un = np.dot(np.dot(QPH, la.inv(np.dot(Phi, QPH))), b) # $u^{(n)}$

        # set small $u^{(n)}$ to zero

        un[abs(un)< thresh]= 0

        # scale $\epsilon$

        if la.norm(un-ul)/la.norm(un)< np.sqrt(epsilon)/eps_thresh:
            epsilon=epsilon/float(eps_scale)
            
        ul = un # last u = this u: $u^{(n-1)}=u^{(n)}$
        
    return un,iter
    

if __name__=="__main__":

    A=np.array([[1,0,1],[1+1j,3,0]])
    b=np.array([[1],[1+1j]])
    x,k=solve(A,b)
    
    print "x: %s"%x
    print "k: %s"%k

