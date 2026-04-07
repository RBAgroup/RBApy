"""Find protein information from UniProt and manual data."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

from itertools import chain
import pandas
import numpy
import re
from rba.prerba.macromolecule import Protein
from rba.prerba.uniprot_data import UniprotData
from rba.prerba.manual_annotation import (
    CuratedCofactors,CuratedLocationMap, 
    CuratedUnknownProteins, CuratedProteins, 
    CuratedGenes)

class ProteinData(object):
    """
    Interface to protein data.

    Attributes
    ----------
    default_stoichiometry : int or float
        Default protein stoichiometry.
    default_location : str
        Default protein location.
    location_map : dict
        Map from UniProt locations to user locations.

    """
    def __init__(self, input_dir,average_gene_id,macromolecule_location_separator):
        """
        Build object from UniProt and manual curation data.
        Parameters
        ----------
        input_dir : str
            Directory with inputs.
        average_gene_id : str
            Location-inspecific prefix of average proteins.
        macromolecule_location_separator : str
            String separating gene-prefix and location identifier in protein IDs.
        """
        self._uniprot = UniprotData(input_dir)
        self._cofactors = CuratedCofactors(input_dir)
        self._user_ids = CuratedUnknownProteins(input_dir)
        self._location_map = CuratedLocationMap(input_dir)
        self._curated_proteins = CuratedProteins(input_dir)
        self._curated_genes = CuratedGenes(input_dir)
        #self._default_location = self._location_map.data.get('Cytoplasm', ('Cytoplasm',None))[0]
        self._default_location = self.query_location_map(key_location='Cytoplasm',value_location='Cytoplasm',comment="OTHER")
        self._average_gene_id = average_gene_id
        self._macromolecule_location_separator = macromolecule_location_separator
        self._average_id = self.average_protein_id(self._default_location)
        self._check_user_identifiers()

    def query_location_map(self,key_location,value_location,comment):
        """"""
        if comment not in ["DoNotAdd"]:
            if self._location_map.data.get(key_location,None) is not None:
                comment_from_map=self._location_map.data.get(key_location,(None,None))[1]
                if comment in ["SBML"]:
                    if comment_from_map in ["FASTA MACHINERY", "GENOME LOCATION", "ANNOTATION PROTEOME","OTHER"]:
                        self._location_map.remove(key_location)
                        self._location_map.append(key_location, value_location, comment)
                elif comment in ["FASTA MACHINERY"]:
                    if comment_from_map in ["GENOME LOCATION", "ANNOTATION PROTEOME","OTHER"]:
                        self._location_map.remove(key_location)
                        self._location_map.append(key_location, value_location, comment)
                elif comment in ["GENOME LOCATION"]:
                    if comment_from_map in ["ANNOTATION PROTEOME","OTHER"]:
                        self._location_map.remove(key_location)
                        self._location_map.append(key_location, value_location, comment)
                elif comment in ["ANNOTATION PROTEOME"]:
                    if comment_from_map in ["OTHER"]:
                        self._location_map.remove(key_location)
                        self._location_map.append(key_location, value_location, comment)
            else:
                self._location_map.append(key_location, value_location, comment)
        return(self._location_map.data.get(key_location,(key_location,None))[0])

    def map_rba_to_sbml_compartment(self,compartment_to_map):
        """Return sbml compartment ID, corresponding to RBA model compartment"""
        return(self._location_map.map_rba_location_to_sbml_location(compartment_to_map=compartment_to_map))

    def _check_user_identifiers(self):
        invalid_identifiers = [i for i in self._user_ids.data.values()
                               if i and not i.startswith(self._average_gene_id)
                               and not self._uniprot.entry(i)]
        if invalid_identifiers:
            print("")
            print('WARNING: {} are invalid gene identifiers. '
                  'Check data provided provided in {}.'
                  .format(', '.join(invalid_identifiers),
                          self._user_ids._raw_data.filename))

    def create_protein_from_gene_id(self, gene_id,protein_id=None,location=None,default_genome_location="",comment_for_location_map="DoNotAdd"):
        """
        Retrieve protein information.

        Parameters
        ----------
        gene_id : str
            Generic location-inspecific ID, without MACROMOLECULE_LOCATION_SEPARATOR... suffix, 
            from sbml or uniprot 
            
        protein_id : str
            Localised protein ID (gene_id + MACROMOLECULE_LOCATION_SEPARATOR... suffix)
            
        location : str
            Compartment where protein is supposed to be located
        Returns
        -------
        Protein
            Basic protein information. None if protein is unknown,
            average protein or 'spontaneous'.
        """
        #if no protein id is provided constructi it from gene id and location
        if not protein_id:
            if not location:
                location=self._default_location
            else:
                protein_id=gene_id+self._macromolecule_location_separator+location
        #use user defined gene id id from unknown protein helper file, if gene is in there. Otherwise use gene id
        user_id = self._user_ids.data.get(gene_id, gene_id)
        if (self._is_spontaneous_id(user_id) or self._is_average_gene_id(user_id)):
            return None
        # get uniprot id from gene id
        uniprot_id = self._uniprot.entry(user_id)
        if uniprot_id: # if uniprot id can be found from gene id
            if gene_id not in self._curated_genes.data:
                self._curated_genes.append(gene_id, default_genome_location)
            # build and return protein
            protein = self.create_protein_from_uniprot_id(uniprot_id,protein_species=protein_id,location=location,origin=self._curated_genes.data.get(gene_id)['LOCATION'],comment=comment_for_location_map)
            return protein
        else:
            # unknown gene identifier
            if gene_id in ['nan']:
                gene_id=numpy.nan
            if gene_id not in self._user_ids.data:
                self._user_ids.append(gene_id, self._average_gene_id)
            return None

    def _is_spontaneous_id(self, id_):
        return id_ == '' or pandas.isnull(id_)

    def _is_average_protein_id(self, id_):
        return(id_.startswith(self._average_gene_id) and id_[int(len(self._average_gene_id)+5):] in [j[0] for j in self._location_map.data.values()])

    def _is_average_gene_id(self,id_):
        return(id_==self._average_gene_id)

    def create_protein_from_uniprot_id(self, uniprot_id,protein_species=None,location=None,origin=None,comment="OTHER"):
        """
        Obtain all necessary information of protein from annotation and manual curation and construct it.
        If not already in manual curation, keys are added.
        """
        if not protein_species:
            protein_species=uniprot_id

        protein = Protein()
        protein.id = protein_species
        protein.origin = origin
        user_location = self._location_map.data.get(location,None)
        if user_location is None:
            protein.location = self._default_location
        else:
            protein.location = self.query_location_map(key_location=location,value_location=location.replace(' ', '_'),comment=comment)
        if origin is not None:
            protein.origin = self.query_location_map(key_location=origin,value_location=origin.replace(' ', '_'),comment="GENOME LOCATION")

        self._fill_with_manual_info(protein, protein_species)
        self._fill_with_uniprot_info(protein, uniprot_id)
        return protein

    def _fill_with_manual_info(self, protein, uniprot_id):
        if protein.origin is None:
            protein.origin = self._curated_genes.data.get(protein.id.split(self._macromolecule_location_separator)[0])['LOCATION']
        protein.cofactors = self._cofactors.data.get(uniprot_id)

    def _fill_with_uniprot_info(self, protein, uniprot_id):
        uniprot_line = self._uniprot.line(uniprot_id)
        #if uniprot_line is not None:
        if protein.location is None:
            protein.location = self._uniprot_location_list(uniprot_line)[0]
        if protein.cofactors is None:
            protein.cofactors = self._uniprot_cofactors(uniprot_line)
        if not protein.sequence:
            protein.sequence = uniprot_line['Sequence']

    def _uniprot_location_list(self, uniprot_line):
        """Return list of compartments from uniprot annotation, mapped via location map. Add entry to location map if necessary."""
        location_list = self._uniprot.find_location_list(uniprot_line)
        if location_list:
            out=[]
            for loc in location_list:
                user_location=self.query_location_map(key_location=loc.replace(' ', '_'),value_location='',comment='ANNOTATION PROTEOME')
                if user_location == '':
                    out.append(loc.replace(' ', '_'))
                else:
                    out.append(user_location)
        else:
            out = [self._default_location]
        return out
    
    def _uniprot_cofactors(self, uniprot_line):
        cofactors, curation_needed = self._uniprot.find_cofactors(uniprot_line)
        if curation_needed:
            if uniprot_line.name not in self._cofactors.data:
                self._cofactors.append(uniprot_line.name, cofactors)
            cofactors = self._cofactors.data.get(uniprot_line.name, [])
        return cofactors

    def _uniprot_subunits(self, uniprot_line,default_value=1):
        subunits = self._uniprot.find_subunits(uniprot_line)
        if not subunits:
            subunits = default_value
        return subunits

    def reference(self, gene_id):
        """
        Retrieve protein reference.

        Returns
        -------
        tuple (str, numeric)
            (id, stoichiometry) that should be used in protein reference.
            If protein is 'spontaneous', None is returned.

        """
        user_id = self._user_ids.data.get(gene_id, gene_id)
        if self._is_spontaneous_id(user_id):
            return None
        if self._is_average_protein_id(user_id):
            return (user_id, 1)
        
        uniprot_id = self._uniprot.entry(gene_id)
        if uniprot_id:
            stoichiometry = self._stoichiometry(uniprot_id)
            return (gene_id, stoichiometry)
        else:
            # unknown gene identifier
            if gene_id in ['nan']:
                gene_id=numpy.nan
            if gene_id not in self._user_ids.data:
                self._user_ids.append(gene_id, self._average_gene_id)
            return (self._average_gene_id, 1)

    def _stoichiometry(self, uniprot_id):
        return self._uniprot_subunits(self._uniprot.line(uniprot_id))

    def average_protein_id(self, compartment):
        """Return identifier of average protein in given compartment."""
        return self._average_gene_id + self._macromolecule_location_separator +  compartment

    def average_composition(self):
        """Return average composition of proteins."""
        return self._uniprot.average_protein_composition()

    def compartments(self):
        """Return list of compartment identifiers."""
        return list(set([j[0] for j in self._location_map.data.values() if j[0] not in ['',None,numpy.nan]]))

    def compartment(self, compartment):
        """Return user identifier associated with compartment."""
        if compartment in [j[0] for j in self._location_map.data.values()]:
            return compartment
        else:
            return self._location_map.data[compartment][0]

    def list_of_mapped_uniprot_compartments(self,gene):
        """Return list of compartments from uniprot annotation, mapped via location map"""
        self._uniprot_location_list(uniprot_line=self._uniprot.line(self._uniprot.entry(self._user_ids.data.get(gene, gene))))
        complist=self._uniprot.find_location_list(uniprot_line=self._uniprot.line(self._uniprot.entry(self._user_ids.data.get(gene, gene))))
        if not complist:
            complist=[self._default_location]
        feasible_location_keys=[i for i in complist if self._location_map.data.get(i.replace(' ', '_'),(None,None))[0] not in ['',None,numpy.nan]]
        out=list(set([self._location_map.data.get(i.replace(' ', '_'),(None,None))[0] for i in feasible_location_keys]))
        return(out)
    
    def annotation_raw_location(self,gene):
        """Get raw location annotation entry from uniprot"""
        return(self._uniprot.line(self._uniprot.entry(self._user_ids.data.get(gene, gene)))['Subcellular location [CC]'])

    def protein_curation_add_entry(self,gene,reaction,column,entry):
        """Add entry to protein curation file"""
        self._curated_proteins._raw_data.data.loc[(self._curated_proteins._raw_data.data['GENE']==gene)&(self._curated_proteins._raw_data.data["REACTION"]==reaction),column]=entry

    def protein_curation_get_entry(self,gene,reaction,column):
        """Get entry from protein curation file"""
        entry_from_file=self._curated_proteins._raw_data.data.loc[(self._curated_proteins._raw_data.data['GENE']==gene)&(self._curated_proteins._raw_data.data["REACTION"]==reaction),column].values[0]
        return(entry_from_file)

    def update_helper_files(self):
        """Update helper files (if needed)."""
        self._location_map.update_file(sort_by=["DATA TYPE","DATA LOCATION"])
        self._cofactors.update_file()
        self._user_ids.update_file()
        self._curated_proteins.update_file()
        self._curated_genes.update_file()
