
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports

# local imports
from rba.prerba.default_data import DefaultParameters


class Enzyme(object):
    DEF_PARAMS = DefaultParameters()

    def __init__(self, reaction, is_transporter):
        self.reaction = reaction
        self.gene_association = []
        self.composition = []
        self.location_from_reaction = [] # location deduced from SBML reactions
        self.compartments_of_metabolites = [] # locations of metabolites inferred from metabolite id suffixes
        self.imported_metabolites = []
        self.is_transporter = is_transporter
        self.initialize_efficiencies()

    def initialize_efficiencies(self):
        if self.is_transporter:
            if self.imported_metabolites:
                self.forward = self.DEF_PARAMS.transport_aggregate_id(
                    self.reaction
                )
            else:
                self.forward = self.DEF_PARAMS.default_transporter_efficiency
            self.backward = self.DEF_PARAMS.default_transporter_efficiency
        else:
            self.forward = self.backward = self.DEF_PARAMS.default_enzyme_efficiency
