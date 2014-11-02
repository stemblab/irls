# Algorithm from <a href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.154.8360">this paper</a>.

# $\LaTeX$ definitions
# brackets: $$\def\lr#1{\left(#1\right)}$$

$blab.irls = (A, b) ->

    # Imports
    sqrt = Math.sqrt
    cpinv = $blab.cpinv # complex SVD blab

    # Helper function
    abs = (z) -> (sqrt(z.x[k]*z.x[k] + z.y[k]*z.y[k]) for k in [0...z.x.length])

    # Algorithm defaults
    p = 0; # norm
    thresh = 0.001 # set coordinates of u to zero if less than this
    max_eps = 1e-8 # accuracy stopping condition
    max_iter = 200 # iteration stopping condition
    eps_scale = 10 # factor to scale epsilon if eps_thresh condition met
    eps_thresh = 100

    # Initialization
    epsilon = 1
    ul = cpinv(A).dot(b) # "last" u estimate (LMS solution)
    m = ul.x.length # .x for real part (issue with .length)

    for k in [1..max_iter]

        break if epsilon<max_eps

        # $q_i=\lr{\lr{u_i^{(n-1)}}^2+\epsilon}^{p/2-1}$ 
        
        
        # $Q_i=1/w_i$

        # Intermediate results
        a = abs ul
        q = nm.pow(a*a + epsilon, 1-p/2)
        Q = complex nm.diag(q), nm.zeros(m, m)
        QAH = Q.dot A.H
        
        # "new" u estimate
        
        # $u^{(n)}=Q_n\Phi^T \lr{\Phi Q_n \Phi^T}^{-1}b$
        
        
        un = (QAH.dot (A.dot QAH).inv()).dot(b)

        # Set small values of un to zero
        for k in [0...m]
            v = [k, 0]
            if un.getBlock(v, v).abs().x[0][0] < thresh
                un.setBlock(v, v, [[0, 0]]) 

        # Scale epsilon
        if (un - ul).norm2() / un.norm2() < sqrt(epsilon)/eps_thresh
            epsilon /= eps_scale

        # new estimate of u becomes last estimate    
        ul = un

    {ul, k}
