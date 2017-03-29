from itertools import izip
import numpy as np
import scipy.stats as st

def bimodality(a):
    mom2 = st.moment(a, 2)
    mom3 = st.moment(a, 3)
    mom4 = st.moment(a, 4)
    return (mom2, mom4, ((mom3**2/mom2**3) + 1.0)/(mom4/mom2**2))

def print_quantity_distributions(analyzer, walkers, label, calculator, calculator_label, do_bimodality=False):
    local_quantities = [calculator(at, calculator_label) for at in walkers]

    if analyzer.comm is None:
        global_quantities = local_quantities
    else:
        global_quantities = analyzer.comm.gather(local_quantities, root=0)

    if analyzer.comm is None or analyzer.comm.rank == 0:
        print "NSAnalyzer:", label, calculator_label, "var", np.var(global_quantities)
        if do_bimodality:
            print "NSAnalyzer:", label, calculator_label, "mu_2 mu_4 Searle_bimodality", bimodality(np.array(global_quantities).flatten())
        global_histo = np.histogram(global_quantities, bins=20)
        bin_width = global_histo[1][1]-global_histo[1][0]
        for (i, (n, r)) in enumerate(izip(global_histo[0], global_histo[1])):
            print "NSAnalyzer:", label, calculator_label, i, r, n/bin_width
