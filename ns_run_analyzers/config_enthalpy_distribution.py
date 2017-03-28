from ns_analyzers_util import print_quantity_distributions

class NSAnalyzer():
    def __init__(self, comm):
        self.comm = comm
        return

    def config_enthalpy(self, at, label):
        try:
            KE = at.get_kinetic_energy()
        except:
            KE = 0.0
        return at.info['ns_energy'] - KE

    def analyze(self, walkers, iter, label):
        print_quantity_distributions(self, walkers, label, self.config_enthalpy, "config_enthalpy")
