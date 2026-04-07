#!/usr/bin/env python3
"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import argparse
import os
import datetime

import rba
from rba.xml._rbaml_version import __rbaml_version__ as rbaml_version

# package imports

def generate_other_macromolecules_file(output_dir):
    rba.xml.RbaMacromolecules().write(os.path.join(output_dir, 'other_macromolecules.xml'))

def generate_custom_constraints_file(output_dir):
    rba.xml.RbaCustomConstraints().write(os.path.join(output_dir, 'custom_constraints.xml'))

def generate_compartments_file(output_dir):
    rba.xml.RbaCompartments().write(os.path.join(output_dir, 'compartments.xml'))

def generate_dummy_model_file_index(output_dir):
    model_file_index=rba.model.ModelFileIndex()
    model_file_index.parameters={'compartments':'compartments.xml', 
                                     'metabolism':'metabolism.xml', 
                                     'enzymes':'enzymes.xml',
                                     'proteins':'proteins.xml',
                                     'rnas':'rnas.xml',
                                     'dna':'dna.xml',
                                     'other_macromolecules':'other_macromolecules.xml',
                                     'processes':'processes.xml',
                                     'targets':'targets.xml',
                                     'parameters':'parameters.xml',
                                     'custom_constraints':'custom_constraints.xml',
                                     'medium':'medium.tsv'}
    model_file_index.write_to_file(file_path=os.path.join(output_dir,'model_file_index.in'))

def add_half_life_attribute_to_macromolecules(model):
    """Add half_life attribute to macromolecules and add respective half_life parameter to model functions (assumed to be constant infinity)"""
    parameter_functions_in_model=[i.id for i in model.parameters.functions]
    for i in model.proteins.macromolecules:
        if i.half_life is None:
            i.half_life="default_protein_half_life"
    if "default_protein_half_life" not in parameter_functions_in_model:
        model.parameters.functions.append(rba.xml.Function("default_protein_half_life", 'constant',{'CONSTANT': 'inf'}))

    for i in model.rnas.macromolecules:
        if i.half_life is None:
            i.half_life="default_rna_half_life"
    if "default_rna_half_life" not in parameter_functions_in_model:
        model.parameters.functions.append(rba.xml.Function("default_rna_half_life", 'constant',{'CONSTANT': 'inf'}))

    for i in model.dna.macromolecules:
        if i.half_life is None:
            i.half_life="default_dna_half_life"
    if "default_dna_half_life" not in parameter_functions_in_model:
        model.parameters.functions.append(rba.xml.Function("default_dna_half_life", 'constant',{'CONSTANT': 'inf'}))

    for i in model.other_macromolecules.macromolecules:
        if i.half_life is None:
            i.half_life="default_half_life"
    if "default_half_life" not in parameter_functions_in_model:
        model.parameters.functions.append(rba.xml.Function("default_half_life", 'constant',{'CONSTANT': 'inf'}))


def add_subunit_id_to_machinery_constituents(machinery_composition_constituents,model):
    """Add comment with subunit info (incrementing number startng from 1) for all species references which refer to macromolecules"""
    subunit_count=0
    for constituent in machinery_composition_constituents:
        if constituent.comment is None:
            subunit_count+=1
            constituent_id=constituent.species
            if model.proteins.macromolecules.get_by_id(constituent_id) is not None:
                constituent.comment="SU_{}".format(subunit_count)
            elif model.rnas.macromolecules.get_by_id(constituent_id) is not None:
                constituent.comment="SU_{}".format(subunit_count)
            elif model.dna.macromolecules.get_by_id(constituent_id) is not None:
                constituent.comment="SU_{}".format(subunit_count)
            elif model.other_macromolecules.macromolecules.get_by_id(constituent_id) is not None:
                constituent.comment="SU_{}".format(subunit_count)

def update_ids_of_machinery_composition_constituent(machinery_composition_constituents,model):
    """Update species IDs in species references in machinery composition: ID --> ID__loc__<location> """
    metabolic_species_in_model=[i.id for i in model.metabolism.species]
    for constituent in machinery_composition_constituents:
        constituent_id=constituent.species
        if "_loc_" in constituent_id:
            continue
        if model.proteins.macromolecules.get_by_id(constituent_id) is not None:
            if constituent_id not in metabolic_species_in_model:
                if constituent_id.startswith('average_protein_'):
                    constituent.species="average_protein_loc_"+str(model.proteins.macromolecules.get_by_id(constituent_id).compartment)
                else:
                    constituent.species=constituent_id+"_loc_"+str(model.proteins.macromolecules.get_by_id(constituent_id).compartment)
        elif model.rnas.macromolecules.get_by_id(constituent_id) is not None:
            if constituent_id not in metabolic_species_in_model:
                constituent.species=constituent_id+"_loc_"+str(model.rnas.macromolecules.get_by_id(constituent_id).compartment)
        elif model.dna.macromolecules.get_by_id(constituent_id) is not None:
            if constituent_id not in metabolic_species_in_model:
                constituent.species=constituent_id+"_loc_"+str(model.dna.macromolecules.get_by_id(constituent_id).compartment)
        elif model.other_macromolecules.macromolecules.get_by_id(constituent_id) is not None:
            if constituent_id not in metabolic_species_in_model:
                constituent.species=constituent_id+"_loc_"+str(model.other_macromolecules.macromolecules.get_by_id(constituent_id).compartment)

def update_ids_of_process_input_species(process,model):
    """Update macromolecule input-species IDs in production and degradation processings"""
    update_ids_of_processing_input_species(processings=process.processings.productions,model=model)
    update_ids_of_processing_input_species(processings=process.processings.degradations,model=model)

def update_ids_of_processing_input_species(processings,model):
    """Update macromolecule species IDs in processing inputs: ID --> ID__loc__<location> """
    metabolic_species_in_model=[i.id for i in model.metabolism.species]
    for processing in processings._elements:
        for input_species in processing.inputs:
            input_id=input_species.species
            if "_loc_" in input_id:
                continue
            if processing.set=="protein":
                if input_id not in metabolic_species_in_model:
                    if input_id.startswith('average_protein_'):
                        input_species.species="average_protein_loc_"+str(model.proteins.macromolecules.get_by_id(input_id).compartment)
                    else:
                        input_species.species=input_id+"_loc_"+str(model.proteins.macromolecules.get_by_id(input_id).compartment)
            elif processing.set=="rna":
                if input_id not in metabolic_species_in_model:
                    input_species.species=input_id+"_loc_"+str(model.rnas.macromolecules.get_by_id(input_id).compartment)
            elif processing.set=="dna":
                if input_id not in metabolic_species_in_model:
                    input_species.species=input_id+"_loc_"+str(model.dna.macromolecules.get_by_id(input_id).compartment)
            elif processing.set=="other_macromolecules":
                if input_id not in metabolic_species_in_model:
                    input_species.species=input_id+"_loc_"+str(model.other_macromolecules.macromolecules.get_by_id(input_id).compartment)

def update_machinery_composition_ids(model,machinery):
    """Update macromolecule species IDs in machinery composition"""
    update_ids_of_machinery_composition_constituent(machinery_composition_constituents=machinery.machinery_composition.reactants,model=model)
    update_ids_of_machinery_composition_constituent(machinery_composition_constituents=machinery.machinery_composition.products,model=model)

def add_machinery_subunit_ids(model,machinery):
    """
    Add comment attribute with subunit id to species reference in machinery composition reactants, which are macromolecules.
    Only for reactants because products are not part of machinery.
    """
    add_subunit_id_to_machinery_constituents(machinery_composition_constituents=machinery.machinery_composition.reactants,model=model)

def update_process_composition_and_input_ids(model):
    """Update IDs of all macromolecules, mentioned in processes.xml"""
    for process in model.processes.processes:
        update_machinery_composition_ids(model=model,machinery=process.machinery) # machinery components
        update_ids_of_process_input_species(process,model) # inputs to processings

def add_enzyme_subunit_ids(model):
    """Add comment attribute with subunit id to species reference in enzyme composition"""
    for enzyme in model.enzymes.enzymes:
        try:
            add_machinery_subunit_ids(model=model,machinery=enzyme)
        except:
            pass

def update_enzyme_composition_ids(model):
    for enzyme in model.enzymes.enzymes:
        try:
            update_machinery_composition_ids(model=model,machinery=enzyme)
        except:
            pass

def update_target_species(target,model):
    """Update target species if they are not metabolic species (ergo macromolecules): ID --> ID__loc__<location> """

    metabolic_species_in_model=[i.id for i in model.metabolism.species]
    if target.species not in metabolic_species_in_model:
        old_species_id=target.species
        if "_loc_" not in old_species_id:
            if model.proteins.macromolecules.get_by_id(target.species) is not None:
                comp=model.proteins.macromolecules.get_by_id(target.species).compartment
                if old_species_id.startswith('average_protein_'):
                    target.species="average_protein_loc_"+str(comp)
                else:
                    target.species=old_species_id+"_loc_"+str(comp)
            elif model.rnas.macromolecules.get_by_id(target.species) is not None:
                comp=model.rnas.macromolecules.get_by_id(target.species).compartment
                target.species=old_species_id+"_loc_"+str(comp)
            elif model.dna.macromolecules.get_by_id(target.species) is not None:
                comp=model.dna.macromolecules.get_by_id(target.species).compartment
                target.species=old_species_id+"_loc_"+str(comp)
            elif model.other_macromolecules.macromolecules.get_by_id(target.species) is not None:
                comp=model.other_macromolecules.macromolecules.get_by_id(target.species).compartment
                target.species=old_species_id+"_loc_"+str(comp)

def update_target_species_ids(model):
    for target_group in model.targets.target_groups:
        for concentration in target_group.concentrations:
            update_target_species(target=concentration,model=model)
        for production_flux in target_group.production_fluxes:
            update_target_species(target=production_flux,model=model)
        for degradation_flux in target_group.degradation_fluxes:
            update_target_species(target=degradation_flux,model=model)

def update_protein_ids(model):
    """GeneID --> GeneID__loc__<location> for all proteins in proteins xml strucutre"""
    metabolic_species_in_model=[i.id for i in model.metabolism.species]
    for i in model.proteins.macromolecules:
        if "_loc_" in i.id:
            continue
        if i.id not in metabolic_species_in_model:
            old_id=i.id
            if old_id.startswith('average_protein_'):
                i.id="average_protein_loc_"+str(i.compartment)
            else:
                i.id=old_id+"_loc_"+str(i.compartment)

def update_rna_ids(model):
    """GeneID --> GeneID__loc__<location> for all rnas in rnas xml strucutre"""
    metabolic_species_in_model=[i.id for i in model.metabolism.species]
    for i in model.rnas.macromolecules:
        if "_loc_" in i.id:
            continue
        if i.id not in metabolic_species_in_model:
            old_id=i.id
            i.id=old_id+"_loc_"+str(i.compartment)

def update_dna_ids(model):
    """DnaID --> DnaID__loc__<location> for all dnas in dna xml strucutre"""
    metabolic_species_in_model=[i.id for i in model.metabolism.species]
    for i in model.dna.macromolecules:
        if "_loc_" in i.id:
            continue
        if i.id not in metabolic_species_in_model:
            old_id=i.id
            i.id=old_id+"_loc_"+str(i.compartment)

def update_other_macromolecule_ids(model):
    """MacromoleculeID --> MacromoleculeID__loc__<location> for all other macromolecules in other macromolecules xml strucutre"""
    metabolic_species_in_model=[i.id for i in model.metabolism.species]
    for i in model.other_macromolecules.macromolecules:
        if "_loc_" in i.id:
            continue
        if i.id not in metabolic_species_in_model:
            old_id=i.id
            i.id=old_id+"_loc_"+str(i.compartment)

def update_macromolecule_ids(model):
    """Update all occurences of all types of macromolecules"""
    update_process_composition_and_input_ids(model=model) #in processes xml structure
    update_enzyme_composition_ids(model=model) #in enzymes xml structure
    update_target_species_ids(model=model) #in targets xml structure
    update_protein_ids(model=model) #in proteins xml structure
    update_rna_ids(model=model) #in rnas xml structure
    update_dna_ids(model=model) #in dna xml structure
    update_other_macromolecule_ids(model=model) #in other macromolecule xml structure

def update_metadata(model,macromolecule_ids_not_updated):
    """Fille model metadata if applicable"""
    if not model.meta_data["Date_first_model_generation"]:
        model.meta_data["RBApy_version_model_generation"]="N/A"
        model.meta_data["Date_first_model_generation"]="N/A"
        model.meta_data["Date_last_model_generation"]="N/A"
        model.meta_data["Date_last_model_generation"]="N/A"
        model.meta_data["Date_download_proteome_annotation"]="N/A"
    model.meta_data["Date_of_model_update_to_v3_or_higher"]=str(datetime.datetime.now())
    model.meta_data["RBApy_version_model_update"]=rba.__version__
    model.meta_data["RBAML_version"]=rbaml_version
    if not macromolecule_ids_not_updated:
        model.meta_data["Localised_Macromolecules"]=True
        model.meta_data["Macromolecule_location_separator"]="_loc_"
    #return(model)

def update_parameters(model):
    update_michaelis_menten_function_with_hill_coefficient(model=model) # add Hill coefficient of 1 to all michaelis menten functions in model
    update_function_references_in_aggregates(model=model) # add exponent (=1) attribute to each function reference in model aggregates

def update_michaelis_menten_function_with_hill_coefficient(model):
    """Add Hill coefficient of 1 to all michaelis menten functions in model"""
    for function in model.parameters.functions:
        if function.type == "michaelisMenten":
            if len([i for i in function.parameters if i.id=="HILL_COEFFICIENT"])==0:
                function.parameters.append(rba.xml.parameters.Parameter("HILL_COEFFICIENT", 1.0))

def update_function_references_in_aggregates(model):
    """Add exponent (=1) attribute to each function reference in model aggregates"""
    for aggregate in model.parameters.aggregates:
        for function_ref in aggregate.function_references:
            function_ref.exponent = 1.0

def populate_compartments_xml(model,output_dir):
    """
    Move information on model compartments from metabolism and density xml file to compartments.xml file. 
    Also add global compartment capacity constraint to custom_constraints.xml
    """
    global_compartment_constraint=rba.xml.custom_constraints.Constraint()
    global_compartment_constraint.id="Global_occupation"
    global_compartment_constraint.action="add"
    global_compartment_constraint.definition.upper_bound=0.0

    global_rhs_param="amino_acid_concentration" #assumed to be parameter in model, dfined in parameters.xml
    global_compartment_constraint.definition.parameter_references.append(rba.xml.custom_constraints.ParameterReference(parameter=global_rhs_param,constant_coefficient=-1.0))

    for compartment in model.metabolism.compartments: 
        corresponding_density=[t for t in model.density.target_densities if t.compartment==compartment.id] 
        if corresponding_density:
            if corresponding_density[0].value is not None:
                model.compartments.compartments.append(rba.xml.compartments.Compartment(id_=compartment.id,lb_=corresponding_density[0].value,ub_=corresponding_density[0].value,is_external_=False))
            elif  corresponding_density[0].upper_bound is not None:
                model.compartments.compartments.append(rba.xml.compartments.Compartment(id_=compartment.id,ub_=corresponding_density[0].upper_bound,is_external_=False))
            elif  corresponding_density[0].lower_bound is not None:
                model.compartments.compartments.append(rba.xml.compartments.Compartment(id_=compartment.id,lb_=corresponding_density[0].lower_bound,is_external_=False))
            else:
                model.compartments.compartments.append(rba.xml.compartments.Compartment(id_=compartment.id,lb_=None,ub_=None,is_external_=False))
            global_compartment_constraint.definition.variable_references.append(rba.xml.custom_constraints.VariableReference(variable="{}_occupation".format(compartment.id),constant_coefficient=1.0))
        else:
            model.compartments.compartments.append(rba.xml.compartments.Compartment(id_=compartment.id,lb_=None,ub_=None,is_external_=True))
    
    model.custom_constraints.constraints.append(global_compartment_constraint)
    
    model.compartments.write(os.path.join(output_dir, 'compartments.xml'))
    model.compartments.write(os.path.join(output_dir, 'custom_constraints.xml'))

    print("")
    print("WARNING: Please review parameter ID, representing total cellular density (e.g. amino-acid concentration) \n " \
    "as parameter-attribute in parameter-reference in Global_occupation constraint in custom_constraints.xml  \n " \
    "and set Global_occupation constraint action-attribute from ignore to add (if constraint is applicable).")
    print("")

def delete_old_model_files(directory):
    old_model_files=['compartments.xml', 
                     'density.xml',
                     'metabolism.xml', 
                     'enzymes.xml',
                     'proteins.xml',
                     'rnas.xml',
                     'dna.xml',
                     'other_macromolecules.xml',
                     'processes.xml',
                     'targets.xml',
                     'parameters.xml',
                     'custom_constraints.xml',
                     'medium.tsv']

    for filename in old_model_files:
        if os.path.isfile(os.path.join(directory, filename)):
            os.remove(os.path.join(directory, filename))

def main():
    """
    Updata rba model from rbaml v1 to rbaml v2
    """
    parser = argparse.ArgumentParser(description='Update a RBA model, to make compatible with RBApy version 3+')
    parser.add_argument('model_dir', metavar='model-dir', type=str,
                        help='Directory for model')
    parser.add_argument('--dont_update_macromolecule_ids', action='store_true',
                        help='Do not add _loc_COMPARTMENT suffix to macromolecule IDs')

    args = parser.parse_args()

    ### Generate empty xml files, necessary for the model to be read by rbapy version>=3 ###
    if not os.path.isfile(os.path.join(args.model_dir, 'other_macromolecules.xml')):
        generate_other_macromolecules_file(output_dir=args.model_dir)

    if not os.path.isfile(os.path.join(args.model_dir, 'custom_constraints.xml')):
        generate_custom_constraints_file(output_dir=args.model_dir)

    compartment_file_generated=False
    if not os.path.isfile(os.path.join(args.model_dir, 'compartments.xml')):
        compartment_file_generated=True
        generate_compartments_file(output_dir=args.model_dir)

    if not os.path.isfile(os.path.join(args.model_dir, 'model_file_index.in')):
        generate_dummy_model_file_index(args.model_dir)

    model = rba.RbaModel.from_xml(args.model_dir) # read model from complete set of xml files #

    # if model already has been updated (info from metadata) abort #
    if model.meta_data["Date_of_model_update_to_v3_or_higher"] is not None:
        print("Model seems to be already updated here: {}".format(model.meta_data["Date_of_model_update_to_v3_or_higher"]))
        return 
    if compartment_file_generated:
        populate_compartments_xml(model=model,output_dir=args.model_dir) # fill compartment xml structure with necessary info

    add_enzyme_subunit_ids(model=model) # add comment with subunit id to enzyme component species refs #

    add_half_life_attribute_to_macromolecules(model=model) # add half life attribute to macromolecules 

    update_parameters(model=model) # update parameter structure

    if not args.dont_update_macromolecule_ids:
        update_macromolecule_ids(model=model) # GeneID --> GeneID__loc__<location>
    
    update_metadata(model=model,macromolecule_ids_not_updated=args.dont_update_macromolecule_ids) # log in metadata of model
    
    model.write(args.model_dir,generate_mean_composition_model=True) # store updated model #
    
    delete_old_model_files(args.model_dir) # delete obsolete model files     

if __name__ == '__main__':
    main()
