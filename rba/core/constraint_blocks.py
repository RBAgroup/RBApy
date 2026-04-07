"""Module defining ConstraintBlocks class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
from rba.core.metabolism import Metabolism
from rba.core.parameters import Parameters
from rba.core.species import Species
from rba.core.enzymes import Enzymes
from rba.core.processes import Processes
from rba.core.targets import Targets
from rba.core.custom_constraints import CustomConstraints
from rba.core.compartments import Compartments


class ConstraintBlocks(object):
    """
    Class converting RBA data into computing substructures.

    Attributes
    ----------
    metabolism : rba.core.metabolism.Metabolism
        Metabolism information.
    parameters : rba.core.parametres.Parameters
        Parameters.
    species : rba.core.species.Species)
        Species information.
    enzymes : rba.core.enzymes.Enzymes
        Enzyme information.
    processes : rba.core.processes.Processes
        Process information.
    targets : rba.core.targets.Targets
        Target information.
    custom_constraints : rba.core.custom_constraints.CustomConstraints
        Custom_constraint information.
    compartments : rba.core.compartments.Compartments
        Compartment information.

    """

    def __init__(self, data):
        """
        Constructor.

        Parameters
        ----------
        data : rba.RbaModel
            RBA model containing raw data.

        """
        # extract metabolism
        self.metabolism = Metabolism(data.metabolism.species,
                                     data.metabolism.reactions)
        # extract parameters
        self.parameters = Parameters(data.parameters.functions,
                                     data.parameters.aggregates)
        # extract base species composition (metabolites + polymers)
        self.species = Species(data, self.metabolism.internal,self.parameters)
        # extract enzyme information
        self.enzymes = Enzymes(data.enzymes, self.species,
                               self.metabolism.reactions, self.parameters)
        # add synthesis reaction for metabolites that are also macromolecules
        (new_reactions, names) = self.species.metabolite_synthesis()
        nb_reactions = len(new_reactions)
        if nb_reactions > 0:
            self.metabolism.add_reactions(new_reactions, names,
                                          [False] * nb_reactions)
        # extract process information
        self.processes = Processes(data.processes.processes,
                                   self.species, self.parameters)
        # extract target information
        self.targets = Targets(data.targets,
                               self.species, self.parameters)
        # extract custom_constraint information
        self.custom_constraints = CustomConstraints(data.custom_constraints,self.parameters)

        # extract compartment information
        self.compartments = Compartments(data.compartments,self.species, self.parameters)

        # setup medium
        self.set_medium(data.medium)

        if self.species.growth_rate_dependent_species_cost:
            print("")
            print("WARNING: Long computation time due to macromolecular species composition is defiend by growth_rate dependent parameters: {}.".format(" , ".join(self.species.growth_rate_dependent_species_cost_parameters)))
            print("Consider removing growth rate as independent variable in parameter definition to speed up computations.")
            
    def set_medium(self, medium):
        """
        Change medium composition.

        Args:
            medium: dict mapping metabolite prefixes with their concentration.
        """

        self.medium=medium
        self.metabolism.set_boundary_fluxes(medium)
        self.parameters.update_medium(medium=medium,growth_rate=0.0)
        self.compute_species_composition_matrices()

    def compute_species_composition_matrices(self):
        ### updates of molecular species processing cost, based on processing input fraction parameters
        self.species.construct_species_matrices(parameters=self.parameters)
        self.enzymes.construct_machinery(species=self.species,parameters=self.parameters)
        self.processes.construct_machinery(species=self.species,parameters=self.parameters)
        self.targets.undetermined_targets.construct_target_matrix(species=self.species,parameters=self.parameters)
        self.targets.determined_targets.construct_target_species(species=self.species,parameters=self.parameters)
        self.compartments.compartment_constituent_definitions(species=self.species,parameters=self.parameters)
