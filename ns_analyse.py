#!/usr/bin/env python

import numpy as np, fileinput, itertools

def read_inputs(args, line_skip=0, line_end=None, interval=1):

    inputs = fileinput.input(files=args)

    (n_walkers, n_cull, n_Extra_DOF) = inputs.readline().split()
    n_walkers = int(n_walkers)
    n_cull = int(n_cull)
    n_Extra_DOF = int(n_Extra_DOF)

    Es=[]
    lines = itertools.islice(inputs, line_skip, line_end, interval)
    for line in lines:
        (n_iter, E) = line.split()
        Es.append(float(E))

    return (n_walkers, n_cull, n_Extra_DOF, np.array(Es))

def calc_log_a(n_Es, n_walkers, n_cull, interval=1):
    # log_a = math.log(float(n_walkers)) - math.log(float(n_walkers+n_cull))
    # From SENS paper PRX v. 4 p 031034 (2014) Eq. 3
    i_range = np.array(range(n_Es*interval))
    i_range_mod_n_cull = np.mod(i_range,n_cull)
    i_range_plus_1_mod_n_cull = np.mod(i_range+1,n_cull)
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
    log_a = log_X_n[0::interval] - np.log(n_walkers+1-i_range_plus_1_mod_n_cull[0::interval])
    return log_a

def calc_Z_terms(beta, log_a, Es):
    #DEBUG for i in range(len(log_a)):
        #DEBUG print "calc_Z_terms log_a ", log_a[i], beta*Es[i]
    shift = np.amax(log_a[:] - beta*Es[:])
    Z_term = np.exp(log_a[:] - beta*Es[:] - shift)
    return (Z_term, shift)
