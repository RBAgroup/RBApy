"""Module defining Density class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy
from scipy.sparse import csr_matrix , coo_matrix , hstack
# local imports
from rba.core.parameter_vector import ParameterVector


class Compartments(object):
    """
    Class computing density-related substructures.

    Attributes
    ----------
    compartments : list of str
        compartment identifiers. compartments[i] is
        the identifier of the compartment involved in the ith constraint.
    compartment_indices : list of int
        compartment_indices[i] is the index (in the list of ids provided at
        construction) of the compartment involved in the ith constraint.
    signs : list of str
        signs of constraints ('E' for equality, 'L' for
        inequality - Lower than)
    values : rba.core.parameter_vector.ParameterVector
        Right-hand side of constraints (depending on μ).

    """

    def __init__(self, compartments, species, parameters):
        """
        Constructor.

        Parameters
        ----------
        compartments : rba.xml.RbaCompartments
            Structure holding compartment information.
        species : rba.core.species.Species
            Species information.
        parameters : rba.core.parameters.Parameters
            Parameter information.

        """
        
        # extract target densities
        self.compartments=[i for i in compartments.compartments if not i.is_external] # only internal compartments
        self.compartment_ids = [i.id for i in self.compartments] # only internal compartments 
        self.compartment_indices = [i[0] for i in enumerate(self.compartment_ids)] # only internal compartments
        self.occupation_rhs = [0.0]*len(self.compartments) # only internal compartments
        self.occupation_signs = ["E"]*len(self.compartments) # only internal compartments
        self.occupation_lb = ParameterVector([i.lower_bound if i.lower_bound else "default_function_CONSTANT_ZERO" for i in self.compartments], parameters)
        self.occupation_ub = ParameterVector([i.upper_bound if i.upper_bound else "default_UB" for i in self.compartments], parameters)

        self.compartment_constituent_fractions(parameters=parameters)
 
        self.number_metabolites=species.production.shape[0]
        self.number_processes=species.prod_proc_cost.shape[0]

        # only relavant if compartment composition is defined (not fully tested)#
        self.compartment_constituent_definitions(species=species,parameters=parameters)


    def compute_CxC_diagonal(self):
        """Construct diagonal matrix with coefficients representing non-housekeeping mass fraction"""
        return(numpy.array([1.0-sum(self.constituent_fractions[compartment].compute()) for compartment in self.compartment_ids]))


    def compute_compartment_composition_metabolite_rows(self,mu):
        """Matrix block corresponding to metabolites used by compartment housekeeping macromolecules"""
        metabolite_compartment_blocks={}
        for comp in self.compartment_ids:
            metabolite_compartment_blocks[comp]=[]
            for i in self.constituent_species[comp]:
                if i["Weight"] !=0:
                    prefactor=i["Fraction"].compute()[0]/i["Weight"]
                    if i["Half_life"] is not None:
                        degradation_rate=numpy.log(2)/i["Half_life"].compute()[0]
                    else:
                        degradation_rate=0.0
                    metabolite_compartment_blocks[comp].append(i["Production"]*prefactor*(mu+degradation_rate))
                    metabolite_compartment_blocks[comp].append(i["Degradation"]*prefactor*degradation_rate)
            if len(metabolite_compartment_blocks[comp])==0:
                metabolite_compartment_blocks[comp].append(csr_matrix((self.number_metabolites, 1)))
        return(coo_matrix(hstack([sum(metabolite_compartment_blocks[comp]) for comp in self.compartment_ids if metabolite_compartment_blocks[comp]])))


    def compute_compartment_composition_process_rows(self,mu):
        """Matrix block corresponding to processing demands used by compartment housekeeping macromolecules"""
        process_compartment_blocks={}
        for comp in self.compartment_ids:
            process_compartment_blocks[comp]=[]
            for i in self.constituent_species[comp]:
                if i["Weight"] !=0:
                    prefactor=i["Fraction"].compute()[0]/i["Weight"]
                    if i["Half_life"] is not None:
                        degradation_rate=numpy.log(2)/i["Half_life"].compute()[0]
                    else:
                        degradation_rate=0.0
                    process_compartment_blocks[comp].append(i["Production_process_cost"]*prefactor*(mu+degradation_rate))
                    process_compartment_blocks[comp].append(i["Degradation_process_cost"]*prefactor*degradation_rate)
            if len(process_compartment_blocks[comp])==0:
                process_compartment_blocks[comp].append(csr_matrix((self.number_processes, 1)))
        return(coo_matrix(hstack([sum(process_compartment_blocks[comp]) for comp in self.compartment_ids if process_compartment_blocks[comp]])))
        

    def compartment_constituent_fractions(self,parameters):
        """Derive parameter vector for housekeeping mass fractions in compartment"""
        self.constituent_fractions={}
        for i in self.compartments:
            constituent_species_fraction=["default_function_CONSTANT_ZERO"]
            for j in i.composition.constituents:
                if j.fraction is not None:
                    constituent_species_fraction.append(j.fraction)
            self.constituent_fractions[i.id]=ParameterVector(constituent_species_fraction, parameters)


    def compartment_constituent_definitions(self,species,parameters):
        """Derive composition and processing demand matrices for compartment housekeeping constituents"""
        self.constituent_species={}
        for i in self.compartments:
            self.constituent_species[i.id]=[]
            for j in i.composition.constituents:
                if j.fraction is not None:
                    constituent_species={}
                    constituent_species["ID"]=j.species
                    constituent_species["Production"]=species.production[:, [species.ids.index(j.species)]].tocsr()
                    constituent_species["Production_process_cost"]=species.prod_proc_cost[:, [species.ids.index(j.species)]].tocsr()
                    constituent_species["Degradation"]=species.degradation[:, [species.ids.index(j.species)]].tocsr()
                    constituent_species["Degradation_process_cost"]=species.deg_proc_cost[:, [species.ids.index(j.species)]].tocsr()
                    constituent_species["Weight"]=species.weight[:, [species.ids.index(j.species)]].sum()
                    try:
                        constituent_species["Half_life"]=ParameterVector([species.half_lifes[species.ids.index(j.species)]], parameters)
                    except:
                        constituent_species["Half_life"]=None
                    constituent_species["Fraction"]=ParameterVector([j.fraction], parameters)
                    self.constituent_species[i.id].append(constituent_species)



        