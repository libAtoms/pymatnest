from itertools import izip
import numpy as np

def print_quantity_distributions(analyzer, walkers, label, calculator, calculator_label):
    local_quantities = [calculator(at, calculator_label) for at in walkers]

    if analyzer.comm is None:
        global_quantities = local_quantities
    else:
        global_quantities = analyzer.comm.gather(local_quantities, root=0)

    if analyzer.comm is None or analyzer.comm.rank == 0:
        print "NSAnalyzer:", label, calculator_label, "var", np.var(global_quantities)
        global_histo = np.histogram(global_quantities, bins=20)
        bin_width = global_histo[1][1]-global_histo[1][0]
        for (i, (n, r)) in enumerate(izip(global_histo[0], global_histo[1])):
            print "NSAnalyzer:", label, calculator_label, i, r, n/bin_width
