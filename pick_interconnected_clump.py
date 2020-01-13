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
def pick_interconnected(rng, N_at, i_list, j_list, N, max_n_attempts=3, print_prefix=""):
    at_list = []

    n_attempts=0
    while len(at_list) < N and n_attempts < max_n_attempts*N:
        n_attempts += 1
        if len(at_list) == 0: # pick a new seed atom
            at_list = [rng.int_uniform(0,N_at)]
            if N == 1:
                return at_list
            ##RNG at_list = [random.randint(0,len(at)-1)]
            i_neighbors = j_list[np.where(i_list == at_list[0])[0]].tolist()
            rng.shuffle_in_place(i_neighbors)
            ##RNG random.shuffle(i_neighbors)

        at_list = add_one_interconnected(i_list, j_list, i_neighbors, at_list)

    if len(at_list) < N:
        return None
    else:
        return at_list
