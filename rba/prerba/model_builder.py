"""Model used to build RBA XML model."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import itertools
from os.path import join
import datetime

# local imports
from rba import RbaModel
from rba.prerba.user_data import UserData
from rba.prerba.enzyme import Enzyme
from rba.prerba.default_targets import DefaultTargets
from rba.prerba.process_builder import ProcessBuilder
from rba.xml._rbaml_version import __rbaml_version__ as rbaml_version

import rba.xml


class ModelBuilder(object):
    """Build a RBA model from user data."""

    def __init__(self, parameter_file, verbose=False):
        """Constructor."""
        self.data = UserData(parameter_file, verbose=verbose)


    def build_model(self):
        """Build and return entire RbaModel."""
        # all unique metabolic reactions in original metabolic network:
        original_metabolic_reactions=set([r.id.split("_duplicate")[0] for r in self.data.sbml_data.reactions])
        model = RbaModel()
        # order matters
        model.enzymes = self.build_enzymes() # test mpved here
        model.metabolism = self.build_metabolism()
        model.processes = self.build_processes() # must be called before building proteins, dnas, rnas and macromolecules
        self.data.protein_data.update_helper_files()
        model.proteins = self.build_proteins()
        model.rnas = self.build_rnas()
        model.dna = self.build_dna()
        model.other_macromolecules = self.build_other_macromolecules()
        model.compartments = self.build_compartments()# must be called after building proteins and processes 
        model.targets = self.build_targets()# must be called after building proteins and processes 
        model.parameters = self.build_parameters()# must be called after building compartments
        model.custom_constraints = self.build_custom_constraints()
        model.medium = self.build_medium()
        model.output_dir = self.data.output_dir()

        # all unique metabolic reactions in rba model after model creation and isoenzyme elimination by protein curation:
        final_metabolic_reactions=set([r.id.split("_duplicate")[0] for r in self.data.sbml_data.reactions])
        
        # load metadata and fill information
        model.get_metadata(join(model.output_dir, 'metadata.tsv'))
        if model.meta_data["Date_first_model_generation"] is None:
            model.meta_data["Date_first_model_generation"]=str(datetime.datetime.now())
        model.meta_data["Date_last_model_generation"]=str(datetime.datetime.now())
        model.meta_data["RBApy_version_model_generation"]=rba.__version__
        model.meta_data["RBAML_version"]=rbaml_version
        model.meta_data["Date_of_model_update_to_v3_or_higher"]='N/A'
        model.meta_data["RBApy_version_model_update"]='N/A'
        model.meta_data["Localised_Macromolecules"]=True
        model.meta_data["Macromolecule_location_separator"]=self.data.protein_data._macromolecule_location_separator
        model.meta_data["Taxon_ID"]=self.data._organism_id()


        if self.data.uniprot_downloaded:
            model.meta_data["Date_download_proteome_annotation"]=str(datetime.datetime.now())
        else:
            model.meta_data["Date_download_proteome_annotation"]="N/A"

        # check if some unique metabolic reactions are removed by eliiminating all isoenzyme during protein curation:
        for rxn in sorted(list(original_metabolic_reactions.difference(final_metabolic_reactions))):
            # If so warn user.
            print("")
            print("WARNING: All isoenzymes for reaction {} removed. Model might be dysfunctional!".format(rxn))
        return model


    def build_metabolism(self):
        """
        Build metabolism part of RBA model.

        Returns
        -------
        rba.xml.RbaMetabolism
            RBA metabolism model in XML format.

        """
        metabolism = rba.xml.RbaMetabolism()

        metabolism.species = self.data.sbml_species()
        metabolism.reactions = self.data.sbml_reactions()
        metabolism.reactions.append(self._atpm_reaction())
        return metabolism


    def _atpm_reaction(self):
        """Generate atp maintenance reaction, turning over atp in cytoplasm"""
        reaction = rba.xml.Reaction(self.data.default.atpm_reaction, False)
        for m in ['ATP', 'H2O']:
            #chose cytosolic isoform 
            id_ = "{}_{}".format(self.data.metabolite_map[m].sbml_id,self.data.generate_mapping_rba_to_sbml_compartments(rba_compartment=self._cytoplasm()))
            if id_:
                reaction.reactants.append(rba.xml.SpeciesReference(id_, 1))
        for m in ['ADP', 'H', 'Pi']:
            #chose cytosolic isoform 
            id_ = "{}_{}".format(self.data.metabolite_map[m].sbml_id,self.data.generate_mapping_rba_to_sbml_compartments(rba_compartment=self._cytoplasm()))
            if id_:
                reaction.products.append(rba.xml.SpeciesReference(id_, 1))
        return reaction


    def build_custom_constraints(self):
        """
        Build custom_constraints part of RBA model.

        Returns
        -------
        rba.xml.RbaCustomConstraints
            RBA model custom_constraints  in XML format.

        """
        constraints=rba.xml.RbaCustomConstraints() # initiate empty constraints structure

        ### define and add global compartment occupation constraint ###
        global_compartment_constraint=rba.xml.custom_constraints.Constraint()
        global_compartment_constraint.id="Global_occupation"
        global_compartment_constraint.action="add"
        global_compartment_constraint.definition.upper_bound=0.0
        global_compartment_constraint.definition.parameter_references.append(rba.xml.custom_constraints.ParameterReference(parameter="amino_acid_concentration",constant_coefficient=-1.0))
        for c_id in self.data.compartments():
            if c_id != self._external():
                global_compartment_constraint.definition.variable_references.append(rba.xml.custom_constraints.VariableReference(variable="{}_occupation".format(c_id),constant_coefficient=1.0))
        constraints.constraints.append(global_compartment_constraint)
        ###
        #Further constraints can be added below ... 

        return(constraints)


    def build_compartments(self):
        """
        Build compartment part of RBA model.

        Returns
        -------
        rba.xml.RbaCompartments
            RBA model compartments  in XML format.

        """
        compartments=rba.xml.RbaCompartments()
        for c_id in self.data.compartments():
            if c_id != self._external():
                compartments.compartments.append(rba.xml.compartments.Compartment(id_=c_id,ub_='{}_density'.format(c_id),is_external_=False))
            else:
                compartments.compartments.append(rba.xml.compartments.Compartment(id_=c_id,is_external_=True))
        return(compartments)


    def build_parameters(self):
        """
        Build parameter part of RBA model.

        Returns
        -------
        rba.xml.RbaParameters
            RBA parameter model in XML format.

        """
        parameters = rba.xml.RbaParameters()
        for fn in self._all_parameter_functions():
            parameters.functions.append(fn)
        for agg in self._all_parameter_aggregates():
            parameters.aggregates.append(agg)
        return parameters


    def _all_parameter_functions(self):
        return (self.data.default.parameters.default_functions()
                + self._protein_functions()
                + self._target_functions()
                + self._efficiency_functions()
                + self._processing_input_functions()
                )


    def _processing_input_functions(self):
        input_fraction_functions=[]
        for input_fraction_parameter in self.process_builder.processing_input_fractions.keys():
            value=self.process_builder.processing_input_fractions[input_fraction_parameter]
            input_fraction_functions.append(self.data.default.parameters.constant_function_from_id_and_value(id_=input_fraction_parameter, constant=value,variable="Temperature"))
        return(input_fraction_functions)


    def _cytoplasm(self):
        return self.data.compartment('Cytoplasm')


    def _external(self):
        return self.data.compartment('Secreted')


    def _other_compartments(self):
        result = self.data.compartments()
        result.remove(self._cytoplasm())
        result.remove(self._external())
        return result


    def _protein_functions(self):
        return(list([self.data.default.parameters.inverse_average_protein_length(self.data.average_protein_length())]+
                     self.data.default.parameters.protein_fraction_functions(compartments=self.data.compartments())+
                     self.data.default.parameters.protein_pg_fraction_functions(compartments=self.data.compartments())
                    )   
                )


    def _target_functions(self):
        return [
            self.data.default.parameters.metabolite_concentration_function(
                target[1], target[2]
            )
            for target in self.data.metabolite_targets()
        ]


    def _efficiency_functions(self):
        fns=[]
        for e in self.data.transport_enzymes():
            fns += self.data.default.parameters.transport_functions(
                e.reaction, e.imported_metabolites
            )
        return fns


    def _all_parameter_aggregates(self):
        return(list(self.data.default.parameters.default_aggregates() + 
                    self._efficiency_aggregates() + 
                    self.data.default.parameters.density_aggregates(compartments=self.data.compartments())+
                    self.data.default.parameters.pg_protein_aggregates(compartments=self.data.compartments())
                    )
                )


    def _append_aggregates(self, parameters, aggs):
        for agg in aggs:
            parameters.aggregates.append(agg)


    def _efficiency_aggregates(self):
        return [self.data.default.parameters.transport_aggregate(e.reaction, e.imported_metabolites) for e in self.data.transport_enzymes() if e.imported_metabolites]


    def build_proteins(self):
        """
        Build protein part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA protein model in XML format.

        """
        builder = MacromoleculeBuilder(macro_molecule_type="protein")
        # components
        for aa in self.data.default.metabolites.aas:
            builder.add_component(aa, '', 'amino_acid', 1)
        for c in self.data.cofactors():
            builder.add_component(c.chebi, c.name, 'cofactor', 0)
        # save ids
        prot_ids = []
        # enzymatic proteins
        for protein in self.data.enzymatic_localised_proteins:
            if protein.id not in prot_ids:
                builder.add_macromolecule(protein.id, protein.location,protein.composition())
                prot_ids.append(protein.id)

        # average proteins
        average_composition = self.data.average_protein_composition()
        for c in self.data.compartments():
            if self.data.protein_data.average_protein_id(c) not in prot_ids:
                builder.add_macromolecule(self.data.protein_data.average_protein_id(c),c, average_composition)
                prot_ids.append(self.data.protein_data.average_protein_id(c))

        # machinery proteins from processes defined and respective fasta files
        for protein in self.process_machinery_components["Proteins"]:
            if protein.id not in prot_ids: 
                loc = protein.location if protein.location else self._cytoplasm()
                builder.add_macromolecule(protein.id, loc, protein.composition())
                prot_ids.append(protein.id)
        return builder.result


    def build_rnas(self):
        """
        Build RNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA RNA model in XML format.

        """
        builder = MacromoleculeBuilder(macro_molecule_type="rna")
        builder.add_component('A', 'Adenosine residue', 'Nucleotide', 2.9036)
        builder.add_component('C', 'Cytosine residue', 'Nucleotide', 2.7017)
        builder.add_component('G', 'Guanine residue', 'Nucleotide', 3.0382)
        builder.add_component('U', 'Uramine residue', 'Nucleotide', 2.7102)
        # user rnas
        rna_ids = []
        for rna in self.data.trnas:
            if rna.id not in rna_ids:
                builder.add_macromolecule(rna.id, rna.location,rna.composition())
                rna_ids.append(rna.id)

        # average RNA
        for mrna in self.data.default_mrnas.values():
            builder.add_macromolecule(mrna.id, mrna.location,self.data.default.metabolites.relative_composition_mrna)

        # machinery rnas
        for rna in self.process_machinery_components["RNAs"]:
            if rna.id not in rna_ids: 
                loc = rna.location if rna.location else self._cytoplasm()
                builder.add_macromolecule(rna.id, loc, rna.composition())
                rna_ids.append(rna.id)
        return builder.result


    def build_dna(self):
        """
        Build DNA part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            RBA DNA model in XML format.

        """
        builder = MacromoleculeBuilder(macro_molecule_type="dna")
        builder.add_component('A', 'Adenosine residue', 'Nucleotide', 0)
        builder.add_component('C', 'Cytosine residue', 'Nucleotide', 0)
        builder.add_component('G', 'Guanine residue', 'Nucleotide', 0)
        builder.add_component('T', 'Thymine residue', 'Nucleotide', 0)
        for dna in self.data.default_dnas.values():
            builder.add_macromolecule(dna.id, dna.location,self.data.default.metabolites.relative_composition_dna)

        return builder.result


    def build_other_macromolecules(self):
        """
        Build other_macromolecules part of RBA model.

        Returns
        -------
        rba.xml.RbaMacromolecules
            (Empty) RBA model other_macromolecules in XML format.
        """
        builder = MacromoleculeBuilder()
        return builder.result


    def build_enzymes(self):
        """
        Build enzyme part of RBA model.

        Returns
        -------
        rba.xml.RbaEnzymes
            RBA enzyme model in XML format.

        """
        enzymes = rba.xml.RbaEnzymes()
        sbml_enzymes=self.data.sbml_enzymes()
        for e in sbml_enzymes["Enzymes"]:
            enzymes.enzymes.append(self._build_enzyme(enzyme=e,subunit_info=sbml_enzymes["Subunits"][e.reaction]))
        # enzyme corresponding to maintenance ATP reaction
        enzymes.enzymes.append(self._build_enzyme(enzyme=Enzyme(self.data.default.atpm_reaction, False)))
        return enzymes


    def _build_enzyme(self, enzyme,subunit_info=None):
        """Build an enzyme"""
        xml_enzyme = rba.xml.Enzyme(enzyme.reaction + '_enzyme',
                                    enzyme.reaction, enzyme.forward,
                                    enzyme.backward)
        if subunit_info:
            for ref in enzyme.composition:
                xml_enzyme.machinery_composition.reactants.append(rba.xml.SpeciesReference(species=ref[0], stoichiometry=ref[1], comment="SU_{}".format(subunit_info.get(ref[0],0))))
        else:
            for ref in enzyme.composition:
                xml_enzyme.machinery_composition.reactants.append(rba.xml.SpeciesReference(species=ref[0], stoichiometry=ref[1]))
        return xml_enzyme


    def build_processes(self):
        """
        Build process part of RBA model.

        Returns
        -------
        rba.xml.RbaProcesses
            RBA process model in XML format.

        """
        # generate processes and processing maps with PRocessBuilder class
        self.process_builder=ProcessBuilder(data=self.data)
        processes , processing_maps , machinery_components = self.process_builder.create_processes_and_processing_maps()

        self.process_machinery_components=machinery_components
        processes_xml = rba.xml.RbaProcesses()
        for p in processes:
            processes_xml.processes.append(p)
        for m in processing_maps:
            processes_xml.processing_maps.append(m)
        self.data=self.process_builder.data
        return(processes_xml)


    def build_targets(self):
        """
        Build target part of RBA model.

        Returns
        -------
        rba.xml.RbaTargets
            RBA targets in XML format.

        """
        targets = rba.xml.RbaTargets()
        for t in self._all_targets():
            if t is not None:
                targets.target_groups.append(t)
        return targets


    def _all_targets(self):
        def_targ = DefaultTargets(self.data)
        return [
            def_targ.translation(),
            def_targ.transcription(),
            def_targ.replication(),
            #def_targ.rna_degradation(),
            def_targ.metabolite_production(),
            def_targ.macrocomponents(self.data.macrocomponents),
            def_targ.maintenance_atp(self.data.default.atpm_reaction)
        ]


    def build_medium(self):
        # !!! we identify metabolites by their prefix !!!
        # M_glc_p and M_glc_e will be seen as the same metabolite M_glc
        return dict.fromkeys(self.data.external_prefixes(),
                             self.data.default.parameters.default_medium_concentration)


    def export_proteins(self, filename):
        self.data.export_proteins(filename)


class MacromoleculeBuilder(object):
    def __init__(self,macro_molecule_type=None):
        if macro_molecule_type=="protein":
            self.result = rba.xml.RbaProteins()
            self.default_half_life_param="default_protein_half_life"
        elif macro_molecule_type=="rna":
            self.result = rba.xml.RbaRNAs()
            self.default_half_life_param="default_rna_half_life"
        elif macro_molecule_type=="dna":
            self.result = rba.xml.RbaDNA()
            self.default_half_life_param="default_dna_half_life"
        else:
            self.result = rba.xml.RbaMacromolecules()
            self.default_half_life_param="default_half_life"
            
    def add_component(self, id_, name, type_, stoichiometry):
        self.result.components.append(
            rba.xml.Component(id_, name, type_, stoichiometry)
        )

    def add_macromolecule(self, id_, location, composition):
        self.result.macromolecules.append(
            rba.xml.Macromolecule(id_=id_, compartment=location, composition=composition, half_life=self.default_half_life_param)
        )
