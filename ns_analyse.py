#!/usr/bin/env python

import numpy as np, fileinput, itertools, sys

def read_inputs(args, line_skip=0, line_end=None, interval=1):

    inputs = fileinput.input(files=args)

    fields = inputs.readline().split()
    if len(fields) == 3:
        (n_walkers, n_cull, n_Extra_DOF) = fields
        flat_V_prior = "False"
        N_atoms = "0"
    elif len(fields) == 5:
        (n_walkers, n_cull, n_Extra_DOF, flat_V_prior, N_atoms) = fields
    else:
        raise("unknown number of fields %d (not 3 or 5) in first line of energies file" % len(fields))
    n_walkers = int(n_walkers)
    n_cull = int(n_cull)
    n_Extra_DOF = int(n_Extra_DOF)
    if flat_V_prior.lower() == "t" or flat_V_prior.lower() == "true":
        flat_V_prior = True
    elif flat_V_prior.lower() == "f" or flat_V_prior.lower() == "false":
        flat_V_prior = False
    else:
        sys.stderr.write("Unknown flat_V_prior string '%s'\n" % flat_V_prior)
        sys.exit(1)
    N_atoms = int(N_atoms)

    Es=[]
    Vs=[]
    lines = itertools.islice(inputs, line_skip, line_end, interval)
    for line in lines:
        fields = line.split()
        if flat_V_prior: # need volumes
            (n_iter, E, V) = fields[0:3]
        else: # no need
            (n_iter, E) = fields[0:2]
        Es.append(float(E))
        if flat_V_prior:
            Vs.append(float(V))

    if len(Vs) == 0:
        Vs = 1.0
    else:
        Vs = np.array(Vs)
    return (n_walkers, n_cull, n_Extra_DOF, flat_V_prior, N_atoms, np.array(Es), Vs)

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

def calc_Z_terms(beta, log_a, Es, flat_V_prior=False, reweight_delta_P=0.0, N_atoms=0, Vs=1.0):
    #DEBUG for i in range(len(log_a)):
        #DEBUG print "calc_Z_terms log_a ", log_a[i], beta*Es[i]
    if not flat_V_prior:
        Vs = 1.0
    shift = np.amax(float(N_atoms)*np.log(Vs) + log_a[:] - beta*(Es[:] + reweight_delta_P*Vs))
    Z_term = np.exp(float(N_atoms)*np.log(Vs) + log_a[:] - beta*(Es[:] + reweight_delta_P*Vs) - shift)
    return (Z_term, shift)
