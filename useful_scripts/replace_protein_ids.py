#!/usr/bin/env python3
"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import argparse
import os
import pandas

import rba

def main():
    """Replace ProteinID prefices according to provided file.
    Inputs: 
        path to rba model to replace IDs in
        path to file specifiying desired replacements of IDs
        (required columns: ORIGINAL (IDs present in model) REPLACEMENT (new ID to replace original one))
    """
    parser = argparse.ArgumentParser(description='Replace protein ids (gene prefix)')
    parser.add_argument('model_dir', metavar='model-dir', type=str,
                        help='path to directory with model')
    parser.add_argument('replacement_map_path', metavar='replacement-map-path', type=str,
                        help='Path to map with gene ids')

    args = parser.parse_args()

    replacement_map=pandas.read_csv(args.replacement_map_path,sep="\t")
    model = rba.RbaModel.from_xml(args.model_dir)

    gene_location_separator=model.meta_data["Macromolecule_location_separator"]

	# loop on uniprot_id to be replaced
	# uniprot_id may appear in proteins, processes, enzymes and targets
    for row in replacement_map.index:
        original_id=replacement_map.loc[row,"ORIGINAL"]
        replacement_id=replacement_map.loc[row,"REPLACEMENT"]

		# replace in proteins
        for protein in model.proteins.macromolecules:
            gene_prefix=protein.id.split(gene_location_separator)[0]
            if gene_prefix==original_id:
                location=protein.id.split(gene_location_separator)[1]
                protein.id=gene_location_separator.join([replacement_id,location])
                
		# replace in process
        for process in model.processes.processes:
            for subunit in process.machinery.machinery_composition.reactants:
                gene_prefix=subunit.species.split(gene_location_separator)[0]
                if gene_prefix==original_id:
                    location=subunit.species.split(gene_location_separator)[1]
                    subunit.species=gene_location_separator.join([replacement_id,location])
            for subunit in process.machinery.machinery_composition.products:
                gene_prefix=subunit.species.split(gene_location_separator)[0]
                if gene_prefix==original_id:
                    location=subunit.species.split(gene_location_separator)[1]
                    subunit.species=gene_location_separator.join([replacement_id,location])
            for production_processing in process.processings.productions:
                for input_species in production_processing.inputs:
                    gene_prefix=input_species.species.split(gene_location_separator)[0]
                    if gene_prefix==original_id:
                        location=input_species.species.split(gene_location_separator)[1]
                        input_species.species=gene_location_separator.join([replacement_id,location])
            for degradation_processing in process.processings.degradations:
                for input_species in degradation_processing.inputs:
                    gene_prefix=input_species.species.split(gene_location_separator)[0]
                    if gene_prefix==original_id:
                        location=input_species.species.split(gene_location_separator)[1]
                        input_species.species=gene_location_separator.join([replacement_id,location])
                        
		# replace in enzymes
        for enzyme in model.enzymes.enzymes:
            for subunit in enzyme.machinery_composition.reactants:
                gene_prefix=subunit.species.split(gene_location_separator)[0]
                if gene_prefix==original_id:
                    location=subunit.species.split(gene_location_separator)[1]
                    subunit.species=gene_location_separator.join([replacement_id,location])
            for subunit in enzyme.machinery_composition.products:
                gene_prefix=subunit.species.split(gene_location_separator)[0]
                if gene_prefix==original_id:
                    location=subunit.species.split(gene_location_separator)[1]
                    subunit.species=gene_location_separator.join([replacement_id,location])
                    
		# replace in targets
        for target_group in model.targets.target_groups:
            for concentration_target in target_group.concentrations:
                gene_prefix=concentration_target.species.split(gene_location_separator)[0]
                if gene_prefix==original_id:
                    location=concentration_target.species.split(gene_location_separator)[1]
                    concentration_target.species=gene_location_separator.join([replacement_id,location])
            for production_target in target_group.production_fluxes:
                gene_prefix=production_target.species.split(gene_location_separator)[0]
                if gene_prefix==original_id:
                    location=production_target.species.split(gene_location_separator)[1]
                    production_target.species=gene_location_separator.join([replacement_id,location])
            for degradation_target in target_group.degradation_fluxes:
                gene_prefix=degradation_target.species.split(gene_location_separator)[0]
                if gene_prefix==original_id:
                    location=degradation_target.species.split(gene_location_separator)[1]
                    degradation_target.species=gene_location_separator.join([replacement_id,location])
                    
    # write RBA model            
    model.write(args.model_dir,generate_mean_composition_model=True)

if __name__ == '__main__':
    main()

