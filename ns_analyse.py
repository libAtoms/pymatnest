#!/usr/bin/env python

import numpy as np, fileinput

def read_inputs(args, skip=0, last_line=-1):

    inputs = fileinput.input(files=args)

    (n_walkers, n_cull, n_Extra_DOF) = inputs.readline().split()
    n_walkers = int(n_walkers)
    n_cull = int(n_cull)
    n_Extra_DOF = int(n_Extra_DOF)


    Es=[]
    i = 0
    for line in inputs:
        if i >= skip:
            (n_iter, E) = line.split()
            Es.append(float(E))
        i += 1
        if last_line > 0 and i >= last_line:
            break

    return (n_walkers, n_cull, n_Extra_DOF, np.array(Es))

def calc_log_a(n_iter, n_walkers, n_cull):
    # log_a = math.log(float(n_walkers)) - math.log(float(n_walkers+n_cull))
    # From SENS paper PRX v. 4 p 031034 (2014) Eq. 3
    i_range = np.array(range(n_iter))
    i_range_mod_n_cull = np.mod(i_range,n_cull)
    # X_n = \prod_{i=0}^n \frac{N-i\%P}{N+1-i\%P}
    # \log X_n = \sum_{i=0}^n \log (N-i\%P) - \log(N+1-i\%P)
    log_X_n_term = np.log(n_walkers-i_range_mod_n_cull) - np.log(n_walkers+1-i_range_mod_n_cull)
    log_X_n = np.cumsum(log_X_n_term)
    # a_n = X_n - X_(n+1)
    # a_n = \prod_{i=0}^n \frac{N-i\%P}{N+1-i\%P} \left(1 - \frac{N-(n+1)\%P}{N+1-(n+1)\%P}\right)
    # \log(a_n) & = \sum_{i=0}^n \left[\log(N-i\%P)-\log(N+1-i\%P)\right] + \log\left(1-\frac{N-(n+1)\%P}{N+1-(n+1)\%P}\right) \\
    #    & = \log(X_n) + \log\left(1-\frac{N-(n+1)\%P}{N+1-(n+1)\%P}\right) \\
    #    & = \log(X_n) + \log\left( \frac{N+1-(n+1)\%P - N + (n+1)\%P}{N+1-(n+1)\%P} \right) \\
    #    & = \log(X_n) + \log \left(  \frac{1}{N+1-(n+1)\%P} \right) \\
    #    & = \log(X_n) -\log \left(  N+1-(n+1)\%P \right)
    log_a = log_X_n - np.log(n_walkers+1-np.mod(i_range+1,n_cull))
    return log_a

def calc_Z_terms(beta, log_a, Es):
    #DEBUG for i in range(len(log_a)):
        #DEBUG print "calc_Z_terms log_a ", log_a[i], beta*Es[i]
    shift = np.amax(log_a[:] - beta*Es[:])
    Z_term = np.exp(log_a[:] - beta*Es[:] - shift)
    return (Z_term, shift)
