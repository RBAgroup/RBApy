"""Interface to manual annotation."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

import os.path
import warnings
from collections import namedtuple
import pandas

from rba.prerba.curation_data import CurationData
from rba.prerba.uniprot_data import Cofactor

Metabolite = namedtuple('Metabolite', 'name sbml_id stoichiometry concentration')


class CuratedData(object):
    def __init__(self, filename, columns):
        self._raw_data = CurationData(filename, columns)
        self.data = {}
        self._warning = ''
        self.file_already_existed=self._raw_data.file_already_existed

    def update_file(self, sort_by=None):
        if self._raw_data.update_file(sort_by=sort_by) and self._warning:
            print("")
            warnings.warn(self._warning, UserWarning)


class CuratedGenes(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'helper_files/gene_locations.tsv')
        super(CuratedGenes, self).__init__(filename, ['GENE ID','GENOME LOCATION'])
        self.data = {r[0]: {'LOCATION':r[1]} for r in self._raw_data.rows()}

        self._warning = (
            'WARNING: Genes with unknown genome locations have been added to '
            'file {}. Execution will continue with default values.'
            .format(filename)
            )

    def append(self,gene, location):
        """Add ambiguous location."""
        self.data[gene] = {'LOCATION':location}
        self._raw_data.add_row((gene,) +(location,))


class CuratedProteins(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'helper_files/protein_curation.tsv')
        super(CuratedProteins, self).__init__(filename, ['GENE',
                                                         'DISABLED',
                                                         'STOICHIOMETRY',
                                                         'REACTION',
                                                         'SUBUNIT',
                                                         'ISOENZYMES', 
                                                         'RBA LOCATION', 
                                                         'STOICHIOMETRY REVIEW CODE',
                                                         'LOCATION REVIEW CODE',
                                                         'ID IN ANNOTATION',
                                                         'LOCATION OF SBML REACTION', 
                                                         'LOCATION FROM ANNOTATION',
                                                         'RAW LOCATION FROM ANNOTATION',  
                                                         'RAW SUBUNIT STRUCTURE FROM ANNOTATION']
                                                         )
        self.data = {r[0]: {'LOCATION':r[6]} for r in self._raw_data.rows()}

    def append(self,gene, reaction,enzyme, location, subunit=None, loc_annot=None, loc_sbml=None, disabled=0, loc_review_code=0,loc_annot_raw=None,stoichiometry=None,stoichiometry_review_code=None,stoich_annot_raw=None,id_annot=None):
        """Add ambiguous location."""
        self.data[gene] = {'LOCATION':location}
        self._raw_data.add_row((gene,)+
                               (disabled,)+
                               (stoichiometry,)+
                               (reaction,)+
                               (subunit,)+
                               (enzyme,)+
                               (location,)+
                               (stoichiometry_review_code,)+
                               (loc_review_code,)+
                               (id_annot,)+
                               (loc_sbml,)+
                               (loc_annot,)+
                               (loc_annot_raw,)+
                               (stoich_annot_raw,))


class CuratedCofactors(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'helper_files/cofactors.tsv')
        super(CuratedCofactors, self).__init__(
            filename,
            ['ENTRY', 'CHEBI', 'NAME',
             'STOICHIOMETRY', 'UNIPROT ANNOTATION']
            )
        self.data = {}
        for row in self._raw_data.rows():
            self._add_to_data(row[0], Cofactor(*row[1:]))
        self._warning = (
            'WARNING: ambiguous UniProt cofactor notes have been added to '
            'file {}. Execution will continue with default values.'
            .format(filename, CurationData.missing_tag)
            )

    def append(self, entry, cofactors):
        """Add ambiguous cofactors."""
        for c in cofactors:
            self._add_to_data(entry, c)
        self._raw_data.add_rows([(entry,) + c for c in cofactors])

    def _add_to_data(self, entry, cofactor):
        """Add cofactor to data only if it has valid information."""
        # do not include cofactors with missing chebi or stoichiometry
        # as well as cofactors with 0 stoichiometry
        sto = cofactor.stoichiometry
        cofactor_list = self.data.setdefault(entry, [])
        if (pandas.notnull(cofactor.chebi)
                and pandas.notnull(sto)
                and float(sto) > 0):
            cofactor_list.append(cofactor)


class CuratedLocationMap(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'helper_files/location_map.tsv')
        super(CuratedLocationMap, self).__init__(
            filename, ['DATA LOCATION', 'RBA LOCATION', 'DATA TYPE']
            )
        self.data = {r[0]: (r[1],r[2]) for r in self._raw_data.rows()}
        # add mandatory compartments (if they are missing)
        self.data.setdefault('Secreted', ('Secreted','DEFAULT'))
        self._warning = (
            'WARNING: UniProt locations with no user-defined '
            'counterpart have been added to {}.'
            .format(filename)
            )

    def append(self, location, default_value, loc_type=None):
        """Add location without user counterpart."""
        self.data[location] = (default_value,loc_type)
        self._raw_data.add_row((location, default_value, loc_type))

    def remove(self, location):
        self.data.pop(location, None)
        self._raw_data.remove_row_by_index(row_index=int(list(self._raw_data.data[self._raw_data.data['DATA LOCATION']==location].index)[0]))

    def map_rba_location_to_sbml_location(self,compartment_to_map):
        d={r[0]: r[1] for r in self._raw_data.rows() if (r[2]=="SBML") and ("," not in r[0])}
        reverted_d={d[i]:i for i in d.keys()}
        mapped_rba_compartments=[]
        for i in d.values():
            if i in mapped_rba_compartments:
                print("")
                print('WARNING: Multiple (original) SBML compartments mapped to same RBA compartment {} in location map'.format(i))
            else:
                mapped_rba_compartments.append(i)
        return(reverted_d.get(compartment_to_map,None))


class CuratedUnknownProteins(CuratedData):
    def __init__(self, input_dir):
        filename = os.path.join(input_dir, 'helper_files/unknown_proteins.tsv')
        super(CuratedUnknownProteins, self).__init__(
            filename, ['GENE ID SBML', 'GENE ID ANNOTATION']
            )
        if self._raw_data.has_missing_information('GENE ID ANNOTATION'):
            raise UserWarning(filename + ': please fill in the'
                              ' GENE ID ANNOTATION column.')
        self.data = {r[0]: r[1] for r in self._raw_data.rows()}
        self._warning = (
            'WARNING: SBML genes not referenced in UniProt have been added to '
            'file {}. Execution will continue with default values.'
            .format(filename)
            )

    def append(self, gene_id, default_value):
        self.data[gene_id] = default_value
        self._raw_data.add_row((gene_id, default_value))


class CuratedMetabolites(CuratedData):
    def __init__(self, input_dir, known_species):
        filename = os.path.join(input_dir, 'helper_files/metabolites.tsv')
        super(CuratedMetabolites, self).__init__(
            filename, ['ID', 'NAME', 'SBML ID', 'STOICHIOMETRY' , 'CONCENTRATION']
            )
        self.data = {}
        invalid_ids = []
        for id_, name, sbml_id, stoich , conc in self._raw_data.rows():
            if pandas.isnull(sbml_id):
                sbml_id = None
            #elif sbml_id not in known_species:
            elif sbml_id not in [i.rsplit("_",1)[0] for i in known_species]:
                invalid_ids.append(id_)
                sbml_id = None
            if pandas.isnull(conc) or conc == '':
                conc = 0
            if pandas.isnull(stoich) or conc == '':
                stoich = 1
            self.data[id_] = Metabolite(name, sbml_id, float(stoich) ,float(conc))
        if invalid_ids:
            print(
                '{}: {} do not map to valid SBML metabolite ids. '
                'These metabolites will be removed from processes and '
                'production targets.'
                .format(filename, ', '.join(invalid_ids))
                )
        self._warning = (
            'WARNING: unidentified key metabolites have been added to file {}.'
            ' These metabolites will be removed from processes and '
            'production targets.'
            .format(filename)
            )

    def append(self, id_, name, sbml_id, stoichiometry=1, concentration=0):
        """Add unrecognized metabolite."""
        self.data[id_] = Metabolite(name, sbml_id, float(stoichiometry), float(concentration))
        self._raw_data.add_row((id_, name, sbml_id, stoichiometry,concentration))


class CuratedMacrocomponents(CuratedData):
    def __init__(self, input_dir, known_species):
        filename = os.path.join(input_dir, 'helper_files/macrocomponents.tsv')
        super(CuratedMacrocomponents, self).__init__(
            filename, ['TARGET_METABOLITE', 'TARGET_CONCENTRATION']
            )
        self.data = {}
        invalid_ids = []
        for met, flux in self._raw_data.rows():
            if met in known_species:
                self.data[met] = float(flux)
            else:
                invalid_ids.append(met)
        if invalid_ids:
            print(
                '{}: {} are not valid SBML metabolite ids. '
                'These entries will be removed from production targets.'
                .format(filename, ', '.join(invalid_ids))
                )
