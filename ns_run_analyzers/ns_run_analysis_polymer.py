import numpy as np
from mpi4py import MPI
from itertools import izip

class NSAnalyzer():
    def __init__(self, comm):
        self.comm = comm
        self.orig_end_to_end_vector = []
        self.config_counter = 0
        return

    def find_ends(self, at):
        return (np.reshape(np.where(at.arrays['angles'] == '_')[0],(-1,2)))

    def end_to_end_distances(self, at):
        dists = []
        if len(self.orig_end_to_end_vector) == 0:
            first_time=True
        else:
            first_time=False
        for ends in self.find_ends(at):
            cur_pos = np.zeros((3))
            for i_at in range(ends[0], ends[1]):
                cur_pos += at.get_distance(i_at, i_at+1, mic=True, vector=True)
            if first_time:
                self.orig_end_to_end_vector.append(cur_pos/np.linalg.norm(cur_pos))
            dists.append(np.linalg.norm(cur_pos))

        return (dists)

    def end_to_end_vector_dots(self, at):
        vec_dots = []
        for ends in self.find_ends(at):
            cur_pos = np.zeros((3))
            for i_at in range(ends[0], ends[1]):
                cur_pos += at.get_distance(i_at, i_at+1, mic=True, vector=True)
            vec_dots.append(np.dot(self.orig_end_to_end_vector[self.config_counter], cur_pos)/np.linalg.norm(cur_pos))
            self.config_counter += 1
            if self.config_counter == len(self.orig_end_to_end_vector):
                self.config_counter = 0
        return (vec_dots)

    def radii_of_gyration(self, at):
        dists = []
        for ends in self.find_ends(at):
            cur_pos = np.zeros((3))
            r_CoM = np.zeros((3))
            vecs = []
            for i_at in range(ends[0], ends[1]):
                cur_pos = at.get_distance(i_at, i_at+1, mic=True, vector=True)
                r_CoM += cur_pos
                vecs.append(cur_pos)
            r_CoM /= float(ends[1]-ends[0]+1)

            r_g = np.zeros((3))
            for cur_pos in vecs:
                r_g += (cur_pos-r_CoM)**2
            dists.append(np.sqrt(r_g[0] + r_g[1] + r_g[2]))

        return (dists)

    def print_quantity_distributions(self, walkers, label, calculator, calculator_label):
        local_dists = [calculator(at) for at in walkers]
        if self.comm is None:
            global_dists = local_dists
        else:
            global_dists = self.comm.gather(local_dists, root=0)

        if self.comm is None or self.comm.rank == 0:
            global_histo = np.histogram(global_dists, bins=20)
            bin_width = global_histo[1][1]-global_histo[1][0]
            for (n, r) in izip(global_histo[0], global_histo[1]):
                print "NSAnalyzer:",label,calculator_label, r, n/bin_width


    def analyze(self, walkers, iter, label):
        if iter < 0: # startup
            self.print_quantity_distributions(walkers, label+" initial", self.end_to_end_distances, "end_to_end_distance")
            self.print_quantity_distributions(walkers, label+" initial", self.end_to_end_vector_dots, "end_to_end_vector_dots")
            self.print_quantity_distributions(walkers, label+" initial", self.radii_of_gyration, "radius_of_gyration")
        elif iter % 10000 == 0: 
            self.print_quantity_distributions(walkers, label+" %d" % iter, self.end_to_end_distances,"end_to_end_distances")
            self.print_quantity_distributions(walkers, label+" %d" % iter, self.end_to_end_vector_dots,"end_to_end_vector_dots")
            self.print_quantity_distributions(walkers, label+" %d" % iter, self.radii_of_gyration, "radius_of_gyration")
