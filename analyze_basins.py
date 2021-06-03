#!/usr/bin/env python

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--kNN', '-n', type=int, help='number of neighbors for connectivity', default=10)
parser.add_argument('--desc_str', '-d', help='init string for SOAP descriptor for similarity',
                    default='soap n_max=8 l_max=8 cutoff={length_scale*2.5} atom_sigma={length_scale*0.2}')
parser.add_argument('--no_desc_average', dest='desc_average', action='store_false', help='do not add "average" keyword to descriptor string')
parser.add_argument('--desc_length_scale', '-l', type=float, help='length scale (first neighbor) for descriptor', required=True)
parser.add_argument('--metric', choices=['dot_product'], help='metric for distance', default='dot_product')
parser.add_argument('--n_walkers', '-K', type=int, help='number of walkers (for config space volume), required for plotting')
parser.add_argument('--n_removed_per_iter', type=int, help='number of walkers removed per iter '
                                                           '(for config space volume', default=1)

parser.add_argument('--min_barrier', type=float, help='ignore basin if barrier separating it from source basin '
                                                      'is smaller than this value', default=-1.0)
parser.add_argument('--min_in_basin', type=int, help='ignore basin if number of samples in basin is less '
                                                     'than this value', default=0)

parser.add_argument('--energy_field', help='info field for relevant energy', default='ns_energy')
parser.add_argument('--iter_field', help='info field for NS iteration', default='iter')
parser.add_argument('--index', '-i', help='ase.io.read "index" for trajectory config reading', default=':')
parser.add_argument('--min_iter', type=int, help='min NS iter to consider')
parser.add_argument('--max_E', type=float, help='max E to consider')

parser.add_argument('--extra_connections', nargs='+', type=int, help='extra connections to add to nn_matrix', default=[])
parser.add_argument('--merge', action='append', nargs=2, type=int, help='merge children by identifier', default=[])

parser.add_argument('--volume_scale', choices=['actual', 'constant', 'linear_in_E', 'linear_in_iter', 'exp_in_E'],
                    help='how to scale volumes when plotting basins', default='exp_in_E')
parser.add_argument('--filled', '-f', action='store_true', help='fill plotted areas')
parser.add_argument('--plot_file', help='file for plot oflandscape chart', default='landscape_chart.pdf')
parser.add_argument('--plot_label_fields', nargs='+', help='extra fields from atoms.info to put in plot labels', default=[])
parser.add_argument('--plot_smooth', type=float, help='smoothing "timescale", in terms of number of iterations, for plotted basin edges', default=1.0)
parser.add_argument('--plot_clean', action='store_true', help='no legend and simpler x axis')
parser.add_argument('--recalculate', '-r', action='store_true', help='do not use saved intermediate results')
parser.add_argument('--update_intermed', '-u', action='store_true', help='update intermediate files if min_iter or max_E changes effective first iter')
parser.add_argument('--write_basin_configs', '-c', action='store_true', help='write pairs of configurations for every basin split')
parser.add_argument('--verbose', '-v', action='store_true', help='verbose output')
parser.add_argument('traj_file', help='trajectory file')
args = parser.parse_args()

assert args.n_removed_per_iter == 1

import sys
import os
import warnings
import json
import pickle
import numpy as np
import ase.io
from quippy.descriptors import Descriptor
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import save_npz, load_npz
from scipy.sparse.csgraph import connected_components
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
import types

from tqdm import tqdm

length_scale = args.desc_length_scale
desc_str = eval("f'"+args.desc_str+"'")
if args.desc_average:
    desc_str += ' average'
sys.stderr.write(f'using desc_str="{desc_str}"\n')

if desc_str.startswith('_INFO_'):
    desc = types.SimpleNamespace(calc = lambda at : { 'data' : [ at.info[desc_str.replace('_INFO_', '').strip()] ] } )
else:
    desc = Descriptor(desc_str)

first_iter_ind = 0
traj_configs = []
reread_intermed = os.path.exists(args.traj_file+'.intermed') and not args.recalculate
if reread_intermed:
    # read intermediate results from file
    warnings.warn('Loading intermediate results')
    with open(args.traj_file+'.intermed') as fin:
        n_lines = int(fin.readline().strip())
        traj_iters = np.loadtxt(fin, max_rows=n_lines)
        traj_energies = np.loadtxt(fin, max_rows=n_lines)
        traj_descs = np.loadtxt(fin, max_rows=n_lines)
    with open(args.traj_file+'.infos', 'rb') as fin:
        traj_infos = pickle.load(fin)
    traj_iters = traj_iters.astype(int).tolist()

    # filter for args-specified min_iter, max_E
    if args.min_iter is not None and args.min_iter < traj_iters[0]:
        raise RuntimeError(f'Requested min_iter={args.min_iter} but first iter found was {traj_iters[0]}')
    if args.max_E is not None and args.max_E > traj_energies[0]:
        warnings.warn(f'Requested max_E={args.max_E} but first energy found was {traj_energies[0]}')

    if args.min_iter is not None or args.max_E is not None:
        for first_iter, t_E in zip(traj_iters, traj_energies):
            if (args.min_iter is None or first_iter >= args.min_iter) and (args.max_E is None or t_E <= args.max_E):
                break
        first_iter_ind = traj_iters.index(first_iter)
        traj_iters = traj_iters[first_iter_ind:]
        traj_energies = traj_energies[first_iter_ind:]
        traj_descs = traj_descs[first_iter_ind:]
        traj_infos = traj_infos[first_iter_ind:]

if args.write_basin_configs or not reread_intermed:
    # need to read configs, perhaps calculate intermed quantities
    sys.stderr.write('Reading configs and optionally computing descriptors\n')
    if not reread_intermed:
        # init arrays for calculation
        traj_energies = []
        traj_iters = []
        traj_infos = []
        traj_descs = []

    for at in tqdm(ase.io.iread(args.traj_file, index=args.index), desc='Read configs and optionally calc desc'):
        config_iter = at.info[args.iter_field]
        config_E = at.info[args.energy_field]
        if (args.min_iter is not None and config_iter < args.min_iter) or (args.max_E is not None and config_E > args.max_E):
            continue
        if args.write_basin_configs:
            # store configs so we can write out minima of each basin
            traj_configs.append(at)

        # might be here only because we need to read configs, in which case don't 
        # recalculate descriptors etc
        if not reread_intermed:
            traj_energies.append(config_E)
            traj_iters.append(int(config_iter))
            traj_infos.append(at.info)
            # mean takes incoherent average, e.g. if it was per center because average=T wasn't used
            d = np.mean(desc.calc(at)['data'], axis=0)
            traj_descs.append(d)

if (args.update_intermed and first_iter_ind > 0) or not reread_intermed:
    traj_energies = np.asarray(traj_energies)
    traj_descs = np.asarray(traj_descs)

    try:
        warnings.warn('Saving intermediate results (iters, energies, descs, infos)')
        with open(args.traj_file+'.intermed', 'w') as fout:
            fout.write(f'{len(traj_iters)}\n')
            np.savetxt(fout, traj_iters)
            np.savetxt(fout, traj_energies)
            np.savetxt(fout, traj_descs)
        with open(args.traj_file+'.infos', 'wb') as fout:
            pickle.dump(traj_infos, fout)
    except:
        warnings.warn(f'Exception while trying to save intermediate results to {args.traj_file+".intermed"}')

if len(traj_descs.shape) == 1:
    # reshape vector into matrix so sklearn functions are happy
    traj_descs = traj_descs.reshape((-1,1))
sys.stderr.write(f'read {len(traj_iters)} {len(traj_energies)} {len(traj_descs)} {len(traj_configs)} configs\n')

if args.recalculate or not os.path.exists(args.traj_file+'.nn_matrix.npz'):
    if args.metric == 'dot_product':
        # default is Minkowski p=2, which is the same as Euclidean, which is the
        # same as dot product if vectors are normalized, which they are for SOAP
        sys.stderr.write('Calculating NNs\n')
        nns = NearestNeighbors(n_neighbors=args.kNN, metric='minkowski', p=2)
        nns.fit(traj_descs)
        sys.stderr.write('Calculating kNN graph\n')
        nn_matrix = nns.kneighbors_graph()
        # k-NN is not symmetric, but adjacency matrix should be
        nn_matrix += nn_matrix.T
    else:
        raise ValueError(f'Unknown metric {args.metric}\n')

    try:
        save_npz(args.traj_file+'.nn_matrix.npz', nn_matrix)
    except:
        warnings.warn(f'Exception while trying to save nn_matrix to {args.traj_file+".nn_matrix.npz"}')
else:
    warnings.warn('Loading nn_matrix results')
    nn_matrix = load_npz(args.traj_file+'.nn_matrix.npz')
    if first_iter_ind > 0:
        nn_matrix = nn_matrix[first_iter_ind:, first_iter_ind:]
        if args.update_intermed:
            try:
                save_npz(args.traj_file+'.nn_matrix.npz', nn_matrix)
            except:
                warnings.warn(f'Exception while trying to save nn_matrix to {args.traj_file+".nn_matrix.npz"}')

for ii, jj in zip(args.extra_connections[0::2], args.extra_connections[1::2]):
    print('adding extra connection', ii, jj)
    nn_matrix[ii, jj] = 1.0

def component_minima_inds(n_components, component_labels, energies, offset=0):
    """find minimum for each connected component
    Parameters
    ----------
    n_components: int
        number of connected components
    component_labels: list(int)
        component label for each configuration in set
    energies: list(float)
        energies for each configuration in set

    Returns
    -------
    min_E_is: list(float), len(min_E_is) == n_components
        array of indices into energies list corresponding to minimum
        energy for each connected component
    """
    min_E_is = []
    for comp_i in range(n_components):
        component_contents = np.where(component_labels == comp_i)[0]
        t_i = np.argmin(energies[component_contents])
        min_E_is.append(component_contents[t_i] + offset)
    return min_E_is


def similarity(metric, d1, d2):
    if args.metric == 'dot_product':
        return sum(d1*d2)
    else:
        raise ValueError(f'Unknown metric {args.metric}\n')


class Basin():
    def __init__(self, identifier, min_E, n_samples, start_iter_i=None, start_E=None):
        self.children = []
        self.children_iters = []
        self.id = identifier
        self.min_E = min_E
        self.n_samples = n_samples

        self.start_iter_i = start_iter_i
        self.start_E = start_E

        # will be filled in later
        # log_weights index is iteration relative to self.start_iter_i
        self.log_weights = []
        self.end_iter_i = None
        self.end_E = None


    def merge_child_id(self, child_id):
        merged_ids = []
        merged_children = []
        # loop over all children (by iter)
        for child in self.children:
            if child_id is None or child_id == child.id:
                # need to merge this child
                found_child = True
                # merge grandchildren
                merged_ids += child.merge_child_id(None)
                merged_children.append(child)
                # add log weights of merged child to this basin
                for abs_iter in range(child.start_iter_i, child.end_iter_i+1):
                    # log(w1 + w2) = log(exp(log_w1) + exp(log_w2)) = log(f exp(low_w1-log_f) + f exp(log_w2-log_f))
                    #    = log(f) + log(exp(log_w1-log_f) + exp(log_w2-log_f))
                    log_scale = self.log_weights[abs_iter - self.start_iter_i]
                    self.log_weights[abs_iter - self.start_iter_i] = log_scale + np.log(np.exp(self.log_weights[abs_iter - self.start_iter_i] - log_scale) + 
                                                                                        np.exp(child.log_weights[abs_iter - child.start_iter_i] - log_scale))
                merged_ids += [child.id]
        if child_id is not None and len(merged_ids) == 0:
            raise RuntimeError(f'Tried to merge {child_id} into {self.id} but no matching child found')

        # remove merged children (but not deeper - those should have been handled by recursive call 
        # in loop) from self.children and children_iter
        for child in merged_children:
            try:
                child_ind = self.children.index(child)
                del self.children[child_ind]
                del self.children_iters[child_ind]
                break
            except ValueError:
                pass

        return merged_ids


    def __str__(self):
        start_E_str = 'None' if self.start_E is None else f'{self.start_E:.3f}'
        end_E_str = 'None' if self.end_E is None else f'{self.end_E:.3f}'
        l = f'id={self.id} min={self.min_E:.3f} started at {self.start_iter_i} {start_E_str} ends at {self.end_iter_i} {end_E_str}'
        l += f' n_samples {self.n_samples} children (id, start_iter) {[(c.id, start_iter) for c, start_iter in zip(self.children, self.children_iters)]}'
        l += f' log_weight=list({len(self.log_weights)})'
        # l += '['
        # for w in self.log_weights:
            # l += f' {w:.3f}'
        # l += ']'
        return l


    def add_child(self, child_basin):
        self.children.append(child_basin)
        self.children_iters.append(child_basin.start_iter_i)


    def child_barrier(self, child):
        return min(child.start_E - child.min_E, child.start_E - self.min_E)


    def get_log_weight(self, iter_i):
        return self.log_weights[iter_i - self.start_iter_i]


    def get_weight(self, iter_i, log_scale=0.0):
        return np.exp(self.log_weights[iter_i - self.start_iter_i] - log_scale)

class Basins():
    def __init__(self):
        self.by_ident = {}


    def add(self, ident, min_E, n_samples, source_ident=None, split_iter_i=0, split_E=None):
        assert ident not in self.by_ident

        new_basin = Basin(ident, min_E, n_samples, split_iter_i, split_E)
        self.by_ident[ident] = new_basin
        if source_ident is not None:
            self.by_ident[source_ident].add_child(new_basin)


    def end(self, ident, end_iter_i, end_E):
        assert self.by_ident[ident].min_E == end_E
        self.by_ident[ident].end_iter_i = end_iter_i
        self.by_ident[ident].end_E = end_E


    def milestones(self):
        m = {}
        for ident, basin in self.by_ident.items():
            if basin.start_iter_i not in m:
                m[basin.start_iter_i] = [[], []]
            m[basin.start_iter_i][0].append(ident)

            if basin.end_iter_i not in m:
                m[basin.end_iter_i] = [[], []]
            m[basin.end_iter_i][1].append(ident)

        return {k: m[k] for k in sorted(m.keys())}


    def __iter__(self):
        by_iter = self.milestones()
        inds = set()
        for iter_i in range(min(by_iter.keys()), max(by_iter.keys())+1):
            if iter_i in by_iter:
                inds |= set(by_iter[iter_i][0])
            yield inds
            if iter_i in by_iter:
                inds -= set(by_iter[iter_i][1])


    def __str__(self):
        ls = ''
        for ident in self.by_ident:
            ls += f'basin minimum index {ident}, {self.by_ident[ident]}\n'
        ls += 'by milestone\n' + str(self.milestones())
        # ls += '\niterable'
        # for iter_i, idents in enumerate(self):
            # ls += f'\n{iter_i} {idents}'

        return ls


    def store_log_weights(self, iter_i, idents, comp_log_weights):
        for ident, c_log_weight in zip(idents, comp_log_weights):
            self.by_ident[ident].log_weights.append(c_log_weight)


    def merge_minimum_criteria(self, min_barrier, min_count):
        last_merge_ids = [None]
        while len(last_merge_ids) > 0:
            last_merge_ids = []

            # loop over existing basins
            for basin in self.by_ident.values():
                # look for basins that existed initially but do not fulfill some criterion
                if basin.start_iter_i == 0 and basin.n_samples < min_count:
                    # merge its children into it
                    for child in basin.children:
                        last_merge_ids += basin.merge_child_id(child.id)
                    last_merge_ids += [basin.id]

                    if len(last_merge_ids) > 0:
                        # delete from self.by_ident any basins that were merged
                        for merge_id in set(last_merge_ids):
                            del self.by_ident[merge_id]
                        # break and restart since we modified self.by_ident
                        break

                # look for children with low barrier or too few samples
                child_merge_ids = []
                for child in basin.children:
                    if basin.child_barrier(child) < min_barrier or child.n_samples < min_count:
                        child_merge_ids.append(child.id)

                # merge any children found
                for merge_id in child_merge_ids:
                    last_merge_ids += basin.merge_child_id(merge_id)

                if len(last_merge_ids) > 0:
                    # delete from self.by_ident any basins that were merged
                    for merge_id in set(last_merge_ids):
                        del self.by_ident[merge_id]
                    # break and restart since we modified self.by_ident
                    break


    def merge(self, id_source, id_merge):
        merged_ids = self.by_ident[id_source].merge_child_id(self.by_ident[id_merge].id)
        for merged_id in merged_ids:
            del self.by_ident[merged_id]


    def get_overall_vol(self, idents, iter_i):
        log_vols = [self.by_ident[ident].get_log_weight(iter_i) for ident in idents]
        log_scale = max(log_vols)
        scaled_overall_vol = np.sum(np.exp(log_vols-log_scale))
        return scaled_overall_vol, log_scale


    def get_rel_weights(self, iter_i, idents, vol_log_scale):
        weights = { ident : self.by_ident[ident].get_weight(iter_i, log_scale=vol_log_scale) for ident in idents }
        tot_weight = sum([w for w in weights.values()])
        for k in weights:
            weights[k] /= tot_weight
        return weights


    def plot(self, ax, traj_energies, vol_scale):
        # should we explicitly pass relevant attributes of args rather than using global?

        init_V = 1
        final_V = 0.025
        # needed for scaled volume of each iter
        dV = init_V - final_V
        traj_init_E = traj_energies[0]
        traj_final_E = traj_energies[-1]
        traj_init_iter = 0
        traj_final_iter = len(traj_energies)-1
        exp_E_f = (traj_final_E - traj_init_E) / np.log(init_V/final_V)
        exp_E_prefactor = np.exp(traj_init_E / exp_E_f)
        cur_width_l = { 'constant'       : lambda iter_i, overall_vol, log_scale : 1.0,
                        'actual'         : lambda iter_i, overall_vol, log_scale : overall_vol*np.exp(log_scale),
                        'exp_in_E'       : lambda iter_i, overall_vol, log_scale : (exp_E_prefactor * np.exp(-(traj_energies[iter_i]/ exp_E_f))),
                        'linear_in_E'    : lambda iter_i, overall_vol, log_scale : (init_V - dV * (traj_energies[iter_i]-traj_init_E)/(traj_final_E-traj_init_E)),
                        'linear_in_iter' : lambda iter_i, overall_vol, log_scale : (init_V - dV * (iter_i-traj_init_iter)/(traj_final_iter-traj_init_iter)) }
        cur_width = cur_width_l[vol_scale]
        left_align_basins = vol_scale == 'constant'

        iters = iter(self)
        idents = next(iters)
        next_left_edge = 0.0
        geom = {}

        # overall volume and current width
        scaled_overall_vol, vol_log_scale = self.get_overall_vol(idents, 0)

        # total width, basin relative widths and corresponding geometry
        width = cur_width(0, scaled_overall_vol, vol_log_scale)
        plot_rel_widths = self.get_rel_weights(0, idents, vol_log_scale)
        basin_order = []
        for ident in idents:
            geom[ident] = (next_left_edge, next_left_edge + width * plot_rel_widths[ident])
            next_left_edge = geom[ident][1]
            basin_order.append(ident)
        prev_geom = { k: tuple(v) for k, v in geom.items() }

        plot_data = {}

        smooth_coeff = 1.0 / args.plot_smooth

        # make sure plot_data exists for every ident
        for ident in self.by_ident:
            plot_data[ident] = ([], [])

        for iter_i, idents in tqdm(enumerate(iters, 1), desc='Plot'):
            scaled_overall_vol, vol_log_scale = self.get_overall_vol(idents, iter_i)
            rel_weights = self.get_rel_weights(iter_i, idents, vol_log_scale)

            width = cur_width(iter_i, scaled_overall_vol, vol_log_scale)
            # loop over existing basins that have new children this iter
            done_idents = set()
            for ident in idents:
                basin = self.by_ident[ident]
                if iter_i in basin.children_iters:
                    iter_child_basins = [b for b, i in zip(basin.children, basin.children_iters) if i == iter_i]
                    superbasin_idents = [b.id for b in iter_child_basins] + [ident]

                    # new overall relative width, smoothed
                    superbasin_overall_rel_width = sum([rel_weights[sub_ident] for sub_ident in superbasin_idents])
                    superbasin_overall_rel_width = (1.0 - smooth_coeff) * plot_rel_widths[ident] + smooth_coeff * superbasin_overall_rel_width

                    # instantaneous relative widths of components of superbasin
                    superbasin_component_rel_widths = self.get_rel_weights(iter_i, superbasin_idents, vol_log_scale)

                    # find edges of superbasin (new smaller basin + new children), make segments
                    superbasin_left_edge =  0.5 * (prev_geom[ident][0] + prev_geom[ident][1]) - 0.5 * width * superbasin_overall_rel_width

                    # update geom and update list of identifiers that were taken care of
                    next_left_edge = superbasin_left_edge
                    cur_basin_order_i = basin_order.index(ident)
                    for sub_basin_id in superbasin_idents:
                        next_left_edge_pos =  next_left_edge
                        plot_rel_widths[sub_basin_id] = superbasin_overall_rel_width * superbasin_component_rel_widths[sub_basin_id]
                        next_right_edge_pos = next_left_edge + width * plot_rel_widths[sub_basin_id]
                        geom[sub_basin_id] = (next_left_edge_pos, next_right_edge_pos)
                        next_left_edge = geom[sub_basin_id][1]
                        done_idents |= set([sub_basin_id])

                        if ident != sub_basin_id:
                            basin_order.insert(cur_basin_order_i, sub_basin_id)
                            cur_basin_order_i += 1

            # loop over existing basins that do not have new children this iter
            for ident in set(idents) - done_idents:
                # just shrink
                ctr = 0.5 * (prev_geom[ident][0] + prev_geom[ident][1])
                plot_rel_widths[ident] = (1.0 - smooth_coeff) * plot_rel_widths[ident] + smooth_coeff * rel_weights[ident]
                left_edge_pos =  ctr - 0.5 * width * plot_rel_widths[ident]
                right_edge_pos = ctr + 0.5 * width * plot_rel_widths[ident]
                geom[ident] = (left_edge_pos, right_edge_pos)

            # close basins that have ended
            for ended_ident in set(geom.keys()) - set(idents):
                del geom[ended_ident]
                ctr = 0.5 * (prev_geom[ended_ident][0] + prev_geom[ended_ident][1])
                dE = traj_energies[iter_i-1] - traj_energies[iter_i]
                plot_data[ended_ident][0].append([(prev_geom[ended_ident][0], traj_energies[iter_i-1]), (ctr, traj_energies[iter_i-1]-0.25*dE)])
                plot_data[ended_ident][1].append([(prev_geom[ended_ident][1], traj_energies[iter_i-1]), (ctr, traj_energies[iter_i-1]-0.25*dE)])

                basin_order.remove(ended_ident)

            if left_align_basins:
                next_left_edge = 0.0
                for ident in basin_order:
                    left_edge_pos = next_left_edge
                    right_edge_pos = next_left_edge + (geom[ident][1] - geom[ident][0])
                    # do not smooth _again_ when doing left aligned
                    # if ident in geom:
                        # left_edge_pos =  (1.0 - smooth_coeff) * geom[indent][0] + smooth_coeff * left_edge_pos
                        # right_edge_pos = (1.0 - smooth_coeff) * geom[indent][1] + smooth_coeff * right_edge_pos
                    geom[ident] = (left_edge_pos, right_edge_pos)
                    next_left_edge = geom[ident][1]

            # create segments to plot
            for ident in idents:
                if ident in prev_geom and ident in geom:
                    plot_data[ident][0].append([(prev_geom[ident][0], traj_energies[iter_i-1]), (geom[ident][0], traj_energies[iter_i])])
                    plot_data[ident][1].append([(prev_geom[ident][1], traj_energies[iter_i-1]), (geom[ident][1], traj_energies[iter_i])])

            # deepcopy geom
            prev_geom = { k: tuple(v) for k, v in geom.items() }

        for ident_i, ident in enumerate(plot_data):
            basin = self.by_ident[ident]
            E_str = 'None' if basin.start_E is None else f'{basin.min_E:.3f}'
            if args.filled:
                # add segments of left edge from top to bottom
                new_data = [ plot_data[ident][0][0][0] ]
                for segment in plot_data[ident][0]:
                    if segment[0] != new_data[-1]:
                        new_data.append(segment[0])
                    new_data.append(segment[1])
                # add segments of right edge from bottom to top
                new_data.append(plot_data[ident][1][-1][1])
                for segment in reversed(plot_data[ident][1]):
                    if segment[1] != new_data[-1]:
                        new_data.append(segment[1])
                    new_data.append(segment[0])
                plot_data[ident] = [new_data]
            else:
                # concatenate segments for left and right edges
                plot_data[ident] = plot_data[ident][0] + plot_data[ident][1]
            label = f'{ident} {traj_iters[ident]} n={self.by_ident[ident].n_samples} Emin={E_str}'
            label += ' ' + ' '.join([field+'='+(f'{traj_infos[ident][field]:.3f}' if isinstance(traj_infos[ident][field], float) else str(traj_infos[ident][field]) ) for field in args.plot_label_fields])
            ax.add_collection(LineCollection(plot_data[ident], colors=f'C{ident_i}', linewidths=1, linestyle='solid',
                                             facecolor=f'C{ident_i}' if args.filled else None,
                                             label=label))

        if args.volume_scale == 'linear':
            ax.set_xlabel('config space volume')
        elif args.volume_scale == 'constant':
            ax.set_xlabel('relative config space volume')
        else:
            if args.plot_clean:
                ax.set_xlabel(f'scaled config space volume')
            else:
                ax.set_xlabel(f'{args.volume_scale} scaled config space volume')
        ax.set_ylabel('E')
        if not args.plot_clean:
            ax.legend(bbox_to_anchor=(1.1, 1), loc='upper left')
        ax.autoscale()


def component_log_weights(n_labels, component_labels, traj_log_weights, iter_i):
    if traj_log_weights is not None:
        cur_traj_log_weights = traj_log_weights[iter_i:]
        log_weights = []
        for comp_i in range(n_labels):
            comp_inds = np.where(component_labels == comp_i)[0]
            log_scale = cur_traj_log_weights[comp_inds[0]]
            # sum(w) = f * sum(w/f) = f * sum(exp(log(w/f))) = f * sum(exp(log_w-log_f))
            # log(sum(w) = log(f) + log(sum(exp(l_w-f)))
            log_weights.append(log_scale + np.log(np.sum(np.exp(cur_traj_log_weights[comp_inds]-log_scale))))
        log_weights = np.asarray(log_weights)
    else:
        log_weights = [None] * n_labels

    return log_weights


if args.n_walkers is not None:
    traj_iters_factors = np.asarray(traj_iters) - traj_iters[0]
    log_a = np.log(args.n_walkers - (args.n_removed_per_iter - 1)) - np.log(args.n_walkers + 1)
    traj_log_weights = log_a*traj_iters_factors
    # normalize to last (smallest) weight = 2 (so log won't be 0)
    traj_log_weights -= (traj_log_weights[-1] - np.log(0.5))
else:
    traj_log_weights = None

if first_iter_ind > 0:
    warnings.warn('Some argument (min_iter or max_E most likely) is leading to skipping some initial iterations, so any saved basins.pckl file is wrong.  RECALCULATING!')
    args.recalculate = True

if args.recalculate or not os.path.exists(args.traj_file+'.basins.pckl'):
    basins = Basins()

    # initialize connected clusters
    n_components, component_labels = connected_components(nn_matrix, directed=False)
    # calculate the index of the minimum E config for each connected component
    m_inds_of_comp = component_minima_inds(n_components, component_labels, traj_energies)
    for c_ind, m_ind in enumerate(m_inds_of_comp):
        basins.add(m_ind, traj_energies[m_ind], sum(component_labels == c_ind))

    c_log_w = component_log_weights(n_components, component_labels, traj_log_weights, 0)
    basins.store_log_weights(0, m_inds_of_comp, c_log_w)

    for E_i, E in tqdm(enumerate(traj_energies[:-1]), desc='Track basins', total=len(traj_energies)-1):
        # loop over decreasing energies, eliminating each one and checking remaining network
        # removing first (row, col) implicitly assumes that configs are sorted by decreasing energy
        nn_matrix = nn_matrix[1:,1:]

        new_n_components, new_component_labels = connected_components(nn_matrix, directed=False)
        # need to translate from rows/columns in current (smaller) matrix to corresponding
        # values in original set, offset is (E_i + 1)
        m_inds_of_new_comp = component_minima_inds(new_n_components, new_component_labels, traj_energies[E_i + 1:], offset=E_i + 1)

        for new_m_ind in set(m_inds_of_new_comp) - set(m_inds_of_comp):
            # found a new minimum
            # find number of cluster it came from (in previous iter, hence using
            #    component_labels rather than new_component_labels)
            # new_m_ind-(E_i+1) is index into new_component_labels array.  Add 1 to 
            #    make it correct for previous iter component_labels, which was 1 longer
            source_cluster = component_labels[1 + new_m_ind - (E_i + 1)]
            # find minimum config index for the source cluster
            source_m_ind = m_inds_of_comp[source_cluster]

            c_ind = m_inds_of_new_comp.index(new_m_ind)
            basins.add(new_m_ind, traj_energies[new_m_ind], sum(new_component_labels == c_ind), source_m_ind, E_i+1, E)

        for missing_m_ind in set(m_inds_of_comp) - set(m_inds_of_new_comp):
            # doesn't exist now, must have ended last iteration
            basins.end(missing_m_ind, E_i, traj_energies[E_i])

        c_log_w = component_log_weights(new_n_components, new_component_labels, traj_log_weights, E_i+1)
        basins.store_log_weights(E_i+1, m_inds_of_new_comp, c_log_w)

        m_inds_of_comp = m_inds_of_new_comp
        component_labels = new_component_labels

    for m_ind in m_inds_of_comp:
        basins.end(m_ind, len(traj_energies)-1, traj_energies[-1])

    if first_iter_ind == 0 or args.update_intermed:
        # only save if this wasn't a restart from intermediate results, or if update_intermed was explicitly specified
        try:
            with open(args.traj_file+'.basins.pckl', 'wb') as fout:
                pickle.dump(basins, fout)
        except:
            warnings.warn(f'Exception while trying to save nn_matrix to {args.traj_file+".basins.pckl"}')
else:
    warnings.warn('Loading basins')
    with open(args.traj_file+'.basins.pckl', 'rb') as fin:
        basins = pickle.load(fin)

basins.merge_minimum_criteria(min_barrier=args.min_barrier, min_count=args.min_in_basin)

if len(args.merge) > 0:
    print('pre merge basins:')
    print(basins)
    for m_ind_source, m_ind_merge in args.merge:
        basins.merge(m_ind_source, m_ind_merge)

print('basins:')
print(basins)

if args.write_basin_configs:
    with open(f'basins.xyz', 'w') as fout:
        for m_ind, basin in basins.by_ident.items():
            traj_configs[m_ind].info['basin_id'] = basin.id
            ase.io.write(fout, traj_configs[m_ind], format='extxyz')

if traj_log_weights is not None:
    import textwrap

    fig = Figure()
    ax = fig.subplots()
    if not args.plot_clean:
        ax.set_title('\n'.join(textwrap.wrap(f'{args}', width=120)))
    basins.plot(ax, traj_energies, args.volume_scale)
    fig.savefig(args.plot_file, bbox_inches='tight')
