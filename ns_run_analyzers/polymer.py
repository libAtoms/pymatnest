import numpy as np
from ns_analyzers_util import print_quantity_distributions

class NSAnalyzer():
    def __init__(self, comm):
        self.comm = comm
        self.orig_end_to_end_vector = []
        self.first_time = True
        return

    def find_ends(self, at):
        return (np.reshape(np.where(at.arrays['angles'] == '_')[0],(-1,2)))

    def end_to_end_distances(self, at, label):
        dists = []
        for ends in self.find_ends(at):
            cur_pos = np.zeros((3))
            for i_at in range(ends[0], ends[1]):
                cur_pos += at.get_distance(i_at, i_at+1, mic=True, vector=True)
            dists.append(np.linalg.norm(cur_pos))

        return (dists)

    def end_to_end_vector_dots(self, at, label):
        vec_dots = []
        for ends in self.find_ends(at):
            cur_pos = np.zeros((3))
            for i_at in range(ends[0], ends[1]):
                cur_pos += at.get_distance(i_at, i_at+1, mic=True, vector=True)
            if self.first_time:
                self.orig_end_to_end_vector.append(cur_pos/np.linalg.norm(cur_pos))
            vec_dots.append(np.dot(self.orig_end_to_end_vector[self.config_counter], cur_pos)/np.linalg.norm(cur_pos))
            self.config_counter += 1
        return (vec_dots)

    def radii_of_gyration(self, at, label):
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


    def analyze(self, walkers, iter, label):
        if iter < 0: # startup
            if label.find("initial_walk") >= 0:
                self.config_counter = 0
                print_quantity_distributions(self, walkers, label, self.end_to_end_vector_dots, "end_to_end_vector_dot")
            self.config_counter = 0
            print_quantity_distributions(self, walkers, label, self.end_to_end_distances, "end_to_end_distance")
            self.config_counter = 0
            print_quantity_distributions(self, walkers, label, self.radii_of_gyration, "radius_of_gyration")
        elif iter % 10000 == 0: 
            self.config_counter = 0
            print_quantity_distributions(self, walkers, label, self.end_to_end_distances,"end_to_end_distances")
            self.config_counter = 0
            print_quantity_distributions(self, walkers, label, self.radii_of_gyration, "radius_of_gyration")

        self.first_time = False
