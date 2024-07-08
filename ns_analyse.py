#!/usr/bin/env python3

from __future__ import print_function
import numpy as np, fileinput, itertools, sys

def read_inputs(args, line_skip=0, line_end=None, interval=1):

    inputs = fileinput.input(files=args)

    fields = inputs.readline().split()
    if len(fields) == 3:
        (n_walkers, n_cull, n_Extra_DOF) = fields
        flat_V_prior = "False"
        N_atoms = "0"
        variable_n_live = "False"
    elif len(fields) == 5:
        (n_walkers, n_cull, n_Extra_DOF, flat_V_prior, N_atoms) = fields
        variable_n_live = "False"
    elif len(fields) == 6:
        (n_walkers, n_cull, n_Extra_DOF, flat_V_prior, N_atoms, variable_n_live) = fields
    else:
        raise ValueError("unknown number of fields %d (not 3 or 5) in first line of energies file" % len(fields))
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
    if variable_n_live.lower() == "t" or variable_n_live.lower() == "true":
        variable_n_live = True
    elif variable_n_live.lower() == "f" or variable_n_live.lower() == "false":
        variable_n_live = False
    else:
        sys.stderr.write("Unknown variable_n_live string '%s'\n" % variable_n_live)
        sys.exit(1)
    N_atoms = int(N_atoms)

    Es=[]
    Vs=[]
    lines = itertools.islice(inputs, line_skip, line_end, interval)
    n_fields = None
    for line in lines:
        fields = line.split()
        if n_fields is not None and n_fields != len(fields):
            sys.stderr.write(f'Mismatch field # prev {n_fields} cur {len(fields)}, skipping\n')
            continue
        if (len(fields) == 3 and not variable_n_live) or (len(fields) == 4 and variable_n_live):
            try:
                E = float(fields[1])
                V = float(fields[2])
            except:
                continue
            if n_fields is None:
                n_fields = len(fields)
            Es.append(E)
            Vs.append(V)
        elif (len(fields) == 2 and not variable_n_live):
            try:
                E = float(fields[1])
            except:
                continue
            if n_fields is None:
                n_fields = len(fields)
            Es.append(E)
        else: # silently skip lines with problems
            sys.stderr.write("WARNING: input line with problem: number of fields not 2 or 3 (without variable_n_live) or 4 (with variable_n_live), or not floats\n")
            continue

    if len(Vs) == 0:
        Vs = None
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

def calc_Z_terms(beta, log_a, Es, flat_V_prior=False, N_atoms=0, Vs=None):
    #DEBUG for i in range(len(log_a)):
        #DEBUG print "calc_Z_terms log_a ", log_a[i], beta*Es[i]
    log_Z_term = log_a[:] - beta*Es[:]
    if flat_V_prior:
        log_Z_term += float(N_atoms)*np.log(Vs[:])
    shift = np.amax(log_Z_term[:])
    Z_term = np.exp(log_Z_term[:] - shift)
    return (Z_term, shift)
