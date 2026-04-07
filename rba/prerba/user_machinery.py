
# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools

# local imports
from rba.prerba.fasta_parser import RbaFastaParser


class UserMachinery(object):
    def __init__(self, filename, protein_data,location_separator,metabolic_species=[]):
        parser = RbaFastaParser(filename,location_separator=location_separator, protein_data=protein_data,metabolic_species=metabolic_species)
        self.proteins = parser.proteins
        self.rnas = parser.rnas

    def protein_ids(self):
        return [p.id for p in self.proteins]

    def rna_ids(self):
        return [r.id for r in self.rnas]

    def composition(self):
        return {m.id: m.stoichiometry
                for m in itertools.chain(self.proteins, self.rnas)}
    
    def has_nonempty_composition(self):
        if len(self.protein_ids())+len(self.rna_ids())>0:
            return(True)
        else:
            return(False)
