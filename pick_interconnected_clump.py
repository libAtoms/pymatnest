import matscipy.neighbours
import numpy as np

# add one more atom connected to all the others (which are all interconnected) in the seed
def add_one_interconnected(i_list, j_list, i_neighbors, seed):
    add_i = -1
    while add_i < 0 and len(i_neighbors) > 0:
        try_add_i = i_neighbors[-1]
        ok = True
        for j in seed:
            if j not in j_list[np.where(i_list == try_add_i)[0]]:
                ok = False
                break
        if ok:
            seed.append(try_add_i)
            return seed
        else:
            ## i_neighbors = i_neighbors[0:-1]
            i_neighbors.pop()
    return []

# find a clump of N atoms from at, all interconnected to each other.
# try max_n_attempts at each cutoff, starting with r_cut, before giving up and increasing cutoff by factor of r_cut_increment (up to max_r_cut)
def pick_interconnected(rng, at, N, r_cut, max_r_cut=0.0, r_cut_increment=1.1, max_n_attempts=10):
    if max_r_cut == 0.0:
        max_r_cut = r_cut
    at_list = []
    while len(at_list) < N and r_cut <= max_r_cut:
        i_list, j_list = matscipy.neighbours.neighbour_list('ij', at, r_cut)

        n_attempts=0
        while len(at_list) < N and n_attempts < max_n_attempts:
            n_attempts += 1
            if len(at_list) == 0: # pick a new seed atom
                at_list = [rng.int_uniform(0,len(at)-1)]
                ##RNG at_list = [random.randint(0,len(at)-1)]
                i_neighbors = j_list[np.where(i_list == at_list[0])[0]].tolist()
                rng.shuffle_in_place(i_neighbors)
                ##RNG random.shuffle(i_neighbors)

            at_list = add_one_interconnected(i_list, j_list, i_neighbors, at_list)

        if len(at_list) < N:
            # start over with longer cutoff
            at_list = []
            r_cut *= r_cut_increment

    if len(at_list) < N:
        return (None, 0.0)
    else:
        return (at_list, r_cut)
