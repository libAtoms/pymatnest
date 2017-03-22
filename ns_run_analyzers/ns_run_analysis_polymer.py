import numpy as np
from mpi4py import MPI
from itertools import izip

class NSAnalyzer():
    def __init__(self, comm):
        self.comm = comm
        return

    def end_to_end_distances(self, at):
        dists = []
        for ends in np.reshape(np.where(at.arrays['angles'] == '_')[0],(-1,2)):
            v = np.zeros((3))
            for i_at in range(ends[0], ends[1]):
                v += at.get_distance(i_at, i_at+1, mic=True, vector=True)
            dists.append(np.linalg.norm(v))

        return dists

    def print_end_to_end_histo(self, walkers, label):
        local_dists = [self.end_to_end_distances(at) for at in walkers]
        if self.comm is None:
            global_dists = local_dists
        else:
            global_dists = self.comm.gather(local_dists, root=0)

        if self.comm is None or self.comm.rank == 0:
            global_histo = np.histogram(global_dists, n_bins=20)
            for (n, r) in izip(global_histo[0], global_histo[1]):
                print "NSAnalyzer:",label," end to end distance ", r, n

    def analyze(self, walkers, iter):
        if iter < 0: # startup
            self.print_end_to_end_histo(walkers, "initial")
        elif iter % 10000 == 0: 
            self.print_end_to_end_histo(walkers, "%d" % iter)
