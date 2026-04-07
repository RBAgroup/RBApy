"""Module defining FastaEntry and FastaParser classes."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

import os.path

# global imports
from collections import namedtuple
from Bio import SeqIO

# local imports
from rba.prerba.macromolecule import Protein, Rna

# class holding fasta entries
FastaEntry = namedtuple('FastaEntry',
                        'id name set_name origin location stoichiometry sequence')


class RbaFastaParser(object):
    """Parse rba-formatted fasta file."""
    def __init__(self, input_file,location_separator, protein_data=None,metabolic_species=[]):
        self._input_file = input_file
        self._protein_data = protein_data
        self.proteins = []
        self.rnas = []
        if os.path.isfile(input_file):
            try:
                for r in SeqIO.parse(input_file, 'fasta'):
                    self._create_molecule(parse_entry(r),metabolic_species,location_separator)
            except IOError:
                print("")
                raise UserWarning('Please provide file ' + input_file)
            except UserWarning as e:
                print("")
                raise UserWarning(input_file + ': ' + str(e))

    def _create_molecule(self, entry,metabolic_species,location_separator):
        if entry.set_name == 'protein':
            self.proteins.append(self._create_protein(entry,location_separator))
        elif entry.set_name == 'rna':
            molecule = Rna()
            self._initialize_molecule(molecule, entry)
            # make sure that that the id is only changed iby addingf the location suffix if macromolecule is no metabolic species
            if molecule.id not in metabolic_species:
                molecule.id = molecule.id + location_separator + molecule.location
            self.rnas.append(molecule)
        else:
            print("")
            raise UserWarning('Unknown molecule type ' + entry.set_name)

    def _create_protein(self, entry,location_separator):
        result = Protein()
        self._initialize_molecule(result, entry)
        result.cofactors = []
        self._protein_data.query_location_map(key_location=result.location,value_location=result.location,comment="FASTA MACHINERY")
        if self._protein_data:
            self._fill_missing_protein_information(result)
        result.id += location_separator + result.location
        return result

    def _fill_missing_protein_information(self, protein):
        try:
            uniprot_protein=self._protein_data.create_protein_from_uniprot_id(protein.id,origin=protein.origin,location=protein.location,comment="DoNotAdd")
            protein.cofactors = uniprot_protein.cofactors
        except:
            if self._has_missing_info(protein):
                print("")
                print(
                    '{}: protein {} has missing information and does not match'
                    ' a known UniProt id. Please fill in all information or '
                    'adapt identifier.'
                )

    def _print_ids_not_found_in_uniprot(self, ids):
        if ids:
            print("")
            print('WARNING ({}): proteins {} could not be retrieved in '
                  'UniProt.'.format(self._filename, ', '.join(ids)))

    def _initialize_molecule(self, molecule, fasta_record):
        molecule.id = fasta_record.id
        molecule.location = fasta_record.location
        molecule.stoichiometry = fasta_record.stoichiometry
        molecule.sequence = fasta_record.sequence
        molecule.origin = fasta_record.origin

    def _has_missing_info(self, molecule):
        return not (molecule.id and molecule.stoichiometry
                    and molecule.sequence)


def parse_entry(record):
    header = record.description
    try:
        [rba, id_, name, set_, orig, loc, sto] = header.split('|')
        sto = float(sto)
    except ValueError:
        invalid_header(header)
    if rba != 'rba':
        invalid_header(header)
    return FastaEntry(id_, name, set_, orig, loc, sto, str(record.seq))


def invalid_header(line):
    """Raise invalid header exception."""
    print("")
    print('Invalid_header\n\t>{}'.format(line))
    print('Expected\n\t>rba|id|name|set_name|origin|location|stoichiometry')
    raise UserWarning('Invalid input file')
