""" Testing the SNAP potential for Si """
import os.path
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms
from ase.build import bulk
from ase.constraints import StrainFilter
from ase.eos import EquationOfState
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from ase.units import kJ
from pymatgen.analysis.wulff import WulffShape
from pymatgen.io.ase import AseAtomsAdaptor
from quippy.potential import Potential

A_KITTEL = 5.43
CALC = Potential(param_filename="../Si_PRX_GAP/gp_iter6_sparse9k.xml")  # "gap"
OPENKIM_MODELS = [
    # "snap"
    "SNAP_ZuoChenLi_2019_Si__MO_869330304805_000",

    # "snap-q"
    "SNAP_ZuoChenLi_2019quadratic_Si__MO_721469752060_000",

    # purja
    "ThreeBodyBondOrder_PPM_PurjaPunMishin_2017_Si__MO_566683736730_000",

    # "elliott"
    "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003",

    # lee-gga
    "SW_LeeHwang_2012GGA_Si__MO_040570764911_001",

    # lee-lda
    "SW_LeeHwang_2012LDA_Si__MO_517338295712_001",

    # balamane
    "SW_BalamaneHauchShi_2017Brittle_Si__MO_381114941873_003",

    # justo
    "EDIP_JustoBazantKaxiras_1998_Si__MO_958932894036_002",

    # stephenson
    "ThreeBodyCluster_SRS_StephensonRadnySmith_1996_Si__MO_604248666067_000",

    # gong
    "ThreeBodyCluster_Gong_Gong_1993_Si__MO_407755720412_000",

    # balamane-1
    "SW_BalamaneHaliciogluTiller_1992_Si__MO_113686039439_005",

    # wang
    "ThreeBodyBondOrder_WR_WangRockett_1991_Si__MO_081872846741_000",

    # mistriotis
    "MFF_MistriotisFlytzanisFarantos_1989_Si__MO_080526771943_001",

    # kaxiras
    "ThreeBodyCluster_KP_KaxirasPandey_1988_Si__MO_072486242437_000",

    # khor
    "ThreeBodyBondOrder_KDS_KhorDasSarma_1988_Si__MO_722489435928_000",

    # biswas
    "ThreeBodyCluster_BH_BiswasHamann_1987_Si__MO_019616213550_000",

    # sw
    "SW_StillingerWeber_1985_Si__MO_405512056662_006"
]
NAME = "gap"
DFT = False


def main():
    f = open(NAME + ".txt", "w")

    # Bulk, diamond-cubic Si
    si_bulk = bulk("Si", crystalstructure="diamond", a=A_KITTEL, cubic=True)
    si_bulk.calc = CALC
    e_si_bulk = si_bulk.get_potential_energy()
    cell = si_bulk.get_cell()

    # Cohesive energy (Kittel = 4.63 eV/atom)
    # Si atom
    si_atom = Atoms("Si")
    si_atom.calc = CALC
    e_si_atom = si_atom.get_potential_energy()
    e_c = e_si_atom - e_si_bulk / 8
    f.write("Calc Ec = {:.2f} eV/atom, Kittel Ec = 4.63 eV/atom\n".format(e_c))

    # Equation of state
    traj = Trajectory("si_bulk-" + NAME + ".traj", "w")
    for x in np.linspace(0.95, 1.05, 5):
        si_bulk.set_cell(cell * x, scale_atoms=True)
        si_bulk.get_potential_energy()
        traj.write(si_bulk)
    configs = read("si_bulk-" + NAME + ".traj@0:5")
    volumes = [x.get_volume() for x in configs]
    energies = [x.get_potential_energy() for x in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    B_gpa = B / kJ * 1.0e24
    f.write("Calc B = {:.2f} GPa".format(B_gpa))
    eos.plot("si_bulk-eos-" + NAME + ".png")

    # Optimized lattice constant
    si_bulk.set_cell(cell)
    sf = StrainFilter(si_bulk)
    opt = BFGS(sf)
    opt.run(0.0005)
    si_bulk.write("si_bulk-" + NAME + ".vasp")

    # Surface energies
    if not os.path.exists("surface-energies-" + NAME + ".csv"):
        paths = glob("si-surfaces/*.cif")
        data = []
        for path in paths:
            facet = path.split("Si-")[1][:-4]
            si_surf = read(path)
            si_surf.calc = CALC
            e_si_surf = si_surf.get_potential_energy()
            n_si_surf = len(si_surf)
            n_si_bulk = len(si_bulk)
            surf_bulk = n_si_surf / n_si_bulk
            sa = np.linalg.det(si_surf.cell[:2, :2])
            e_s = (e_si_surf - e_si_bulk * surf_bulk) / (2 * sa)
            data.append([facet, e_s])
        df = pd.DataFrame(data, columns=["facet", "e_s"])
        df.to_csv("surface-energies-" + NAME + ".csv", index=False)

    # Wulff shape
    if DFT:
        df = pd.read_excel("surface-energies.xlsx")
    else:
        df = pd.read_csv("surface-energies-" + NAME + ".csv")
    lattice = AseAtomsAdaptor.get_structure(si_bulk).lattice
    miller_list = []
    for x in df.facet:
        y = x[1:4]
        z = (int(y[0]), int(y[1]), int(y[2]))
        miller_list.append(z)
    if DFT:
        e_surf_list_dft = df.dft.values.tolist()
    e_surf_list_gap = df.e_s.values.tolist()

    # DFT Wulff shape
    if DFT:
        wulffshape_dft = WulffShape(lattice, miller_list, e_surf_list_dft)
        plt.figure(figsize=(3.33, 3.33))
        wulffshape_dft.get_plot()
        plt.savefig("dft-wulff-shape.png")

    # GAP
    wulffshape_gap = WulffShape(lattice, miller_list, e_surf_list_gap)
    plt.figure(figsize=(3.33, 3.33))
    wulffshape_gap.get_plot()
    plt.savefig(NAME + "-wulff-shape.png")
    f.close()


if __name__ == "__main__":
    main()
