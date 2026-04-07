"""Module defining ConstraintMatrix class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy
from scipy.sparse import coo_matrix, diags, hstack, vstack

# local imports
from rba.core.constraint_blocks import ConstraintBlocks


class ConstraintMatrix(object):
    """
    Class building constraint matrix.

    Attributes:
        col_names: Linear problem column names (decision variables).
        reaction_cols: Indices of columns corresponding to reactions.
        enzyme_cols: Indices of columns corresponding to enzymes.
        process_cols: Indices of columns corresponding to processes.
        target_cols: Indices of columns corresponding to targets.
        row_names: Linear problem row names (constraints).
        row_signs: Linear problem row signs (equality or inequality).
        UB: Linear problem upper bounds.
        LB: Linear problem lower bounds.
        f: Linear problem objective function.
        A: Linear problem matrix (left-hand side).
        b: Linear problem right-hand side.

    """

    def __init__(self, model):
        """
        Build constraint matrix from model.

        Parameters
        ----------
        model : rba.RbaModel
            RBA model.

        """
        self._blocks = ConstraintBlocks(model)
        # convenience variables
        reactions = self._blocks.metabolism.reactions
        enzymes = self._blocks.enzymes.ids
        processes = self._blocks.processes.ids
        undetermined_fluxes = self._blocks.targets.undetermined_targets.names
        compartments = self._blocks.compartments.compartment_ids
        nb_reactions = len(reactions)
        nb_enzymes = len(enzymes)
        nb_processes = len(processes)
        nb_undetermined = len(undetermined_fluxes)
        nb_compartments = len(compartments)

        # column information
        self.col_names = (reactions                                           # decision variable representing metabolic fluxes
                          + [e for e in enzymes]                              # decision variables representing abundances of enzymes
                          + [p + '_machinery' for p in processes]             # decision variables representing abundances of process-machineries
                          + [m + '_target_flux' for m in undetermined_fluxes]
                          + [c + '_occupation' for c in compartments]) 
        self.reaction_cols = numpy.arange(nb_reactions)
        self.enzyme_cols = nb_reactions + numpy.arange(nb_enzymes)
        self.process_cols = (nb_reactions + nb_enzymes +
                             numpy.arange(nb_processes))
        self.target_cols = (nb_reactions + nb_enzymes + nb_processes +
                            numpy.arange(nb_undetermined))
        # row information
        self.row_names = (self._blocks.metabolism.internal             # Constraints representing metabolite mass-balances
                          + [p + '_capacity' for p in processes]       # Constraints representing occupation and capacity of macromolecular processes
                          + [e + '_forward_capacity' for e in enzymes] # Constraints representing occupation and capacity of enzymes in forward direction
                          + [e + '_backward_capacity' for e in enzymes]# Constraints representing occupation and capacity of enzymes in backward direction
                          + [c + '_density' for c in compartments]     # Constraints representing occupation and density limit of model compartments
                          )       

        self._test_validity_of_custom_constraints()                
        self.row_names += self._blocks.custom_constraints.ids # Custom constraint IDs
        
        self.original_row_names=self.row_names
        self.row_signs = (['E'] * len(self._blocks.metabolism.internal) 
                          + self._blocks.processes.capacity_signs
                          + ['L'] * 2 * nb_enzymes
                          + self._blocks.compartments.occupation_signs
                          + self._blocks.custom_constraints.signs)
        # constant building blocks:
        self._empty_ExPU = coo_matrix((nb_enzymes,
                                       nb_processes + nb_undetermined))
        self._empty_PxR = coo_matrix((nb_processes, nb_reactions))
        self._empty_CxR = coo_matrix((nb_compartments, nb_reactions))
        self._empty_2E = numpy.zeros(2 * nb_enzymes)
        self._empty_ExC = coo_matrix((nb_enzymes,nb_compartments))

        self._eye_CxC = coo_matrix(numpy.eye(nb_compartments))

        # indicator matrices:
        R_ind = [reactions.index(r) for r in self._blocks.enzymes.reaction_catalyzed]
        self._R_to_E = coo_matrix(([1]*nb_enzymes, (range(nb_enzymes), R_ind)),shape=(nb_enzymes, nb_reactions))
        target_reactions = self._blocks.targets.target_reactions
        self._value_reaction_cols = self.reaction_cols[[reactions.index(r) for r in target_reactions.value_reactions]]
        self._lb_reaction_cols = self.reaction_cols[[reactions.index(r) for r in target_reactions.lb_reactions]]
        self._ub_reaction_cols = self.reaction_cols[[reactions.index(r) for r in target_reactions.ub_reactions]]
        # set remaining attributes to None
        self.A = self.b = self.LB = self.UB = self.f = None

    def build_matrices(self, mu):
        """
        Build LP matrices corresponding to given growth-rate.

        Args:
            mu: growth_rate
        """
        # update parameters
        self._blocks.parameters.update_growth_rate(medium=self._blocks.medium, growth_rate=mu)

        if self._blocks.species.growth_rate_dependent_species_cost:
            self._blocks.compute_species_composition_matrices()

        # build A
        enzymes = self._blocks.enzymes
        processes = self._blocks.processes
        targets = self._blocks.targets

        # mu-dependent blocks
        u_composition, u_proc_cost, u_weight = targets.undetermined_targets.matrices(mu)
        process_capacity = processes.capacity.compute()

        forward, backward = enzymes.efficiencies()
        
        # metabolite mass balance constraints:
        metab_rows = _build_metabolite_constraints(stoichiomerty_matrix_metabolism=self._blocks.metabolism.S,
                                 enzymes=enzymes,
                                 processes=processes,
                                 targets_metabolite_components=u_composition,
                                 compartments=self._blocks.compartments,
                                 growth_rate=mu)

        # process-capacity constraints:
        process_rows =  _build_process_capacity_constraints(empty_process_reaction_block=self._empty_PxR,
                                       enzymes=enzymes,
                                       processes=processes,
                                       targets_machinery_cost=u_proc_cost,
                                       process_capacities=process_capacity,
                                       compartments=self._blocks.compartments,
                                       growth_rate=mu)

        # enzyme-capacity constraints:
        forward_rows = hstack([self._R_to_E, -diags(forward), self._empty_ExPU, self._empty_ExC],)
        backward_rows = hstack([-self._R_to_E, -diags(backward), self._empty_ExPU, self._empty_ExC])
        
        # density constraints:
        density_rows = hstack([self._empty_CxR,
                               enzymes.machinery.weight[self._blocks.compartments.compartment_indices],
                               processes.machinery.weight[self._blocks.compartments.compartment_indices],
                               u_weight[self._blocks.compartments.compartment_indices],
                               -coo_matrix(numpy.diag(self._blocks.compartments.compute_CxC_diagonal()))
                               ])
        
        # Stack growth-rate specific metabolite mass-balance, process-capacity, 
        # enzyme forward and backward capacity, compartment density and custom constraints 
        # on top of each other to obtain full µ-specific RBA constraint matrix (LHS):
        self.A = vstack([metab_rows, process_rows,
                         forward_rows, backward_rows, density_rows,self._blocks.custom_constraints.build_custom_constraints_lhs(col_names=self.col_names)])

        # build b (RHS):

        # gather mu-dependent blocks:
        fluxes, processing, weight = targets.determined_targets.compute(mu) 
        #density_rows = density.values.compute() - weight[c_indices].T

        # build vector
        self.b = numpy.concatenate([-fluxes, 
                                    -processing,
                                    self._empty_2E, 
                                    -weight[self._blocks.compartments.compartment_indices],
                                    self._blocks.custom_constraints.build_custom_constraints_rhs()]) 
        
        # update lower bounds and upper bounds:
        self.LB = numpy.concatenate([self._blocks.metabolism.lb(),
                                     self._blocks.enzymes.lb,
                                     processes.lb,
                                     targets.undetermined_targets.lb(),
                                     self._blocks.compartments.occupation_lb.compute()])
        self.UB = numpy.concatenate([self._blocks.metabolism.ub(),
                                     self._blocks.enzymes.ub,
                                     processes.ub,
                                     targets.undetermined_targets.ub(),
                                     self._blocks.compartments.occupation_ub.compute()])
        
        # define objective function
        self.f = numpy.concatenate([self._blocks.metabolism.f,
                                    self._blocks.enzymes.f,
                                    processes.f,
                                    targets.undetermined_targets.f,
                                    numpy.zeros(len(self._blocks.compartments.compartment_ids))])
        
        # target reactions:
        self.LB[self._lb_reaction_cols] = targets.target_reactions.lb()
        self.UB[self._ub_reaction_cols] = targets.target_reactions.ub()
        r_fluxes = targets.target_reactions.value()
        self.LB[self._value_reaction_cols] = r_fluxes
        self.UB[self._value_reaction_cols] = r_fluxes

        # if there are constraints to be removed
        if self._blocks.custom_constraints.constraints_to_remove:
            # remove constraints to be removed
            self.remove_constraints()

    def set_medium(self, medium):
        """
        Change medium composition.

        Args:
            medium: dict mapping metabolite prefixes with their concentration.
        """
        self._blocks.set_medium(medium)

    def remove_constraints(self):
        """
        Removes rows corresponding to constraints to be removed, qs defined in custom constraints from RBA problem.

        """
        row_indices_to_remove=[self.original_row_names.index(i) for i in self._blocks.custom_constraints.constraints_to_remove]
        A_intermediate=self.A.tocsr()
        mask=numpy.ones(self.A.shape[0],dtype=bool)
        mask[numpy.array(row_indices_to_remove)]=False
        self.A=A_intermediate[mask,:].tocoo()

        b_intermediate=[entry for i, entry in enumerate(list(self.b)) if i not in row_indices_to_remove]
        self.b=numpy.array(b_intermediate)

        if len(self.original_row_names)==len(self.row_names):
            row_signs_intermediate=[entry for i, entry in enumerate(self.row_signs) if i not in row_indices_to_remove]
            self.row_signs=row_signs_intermediate
            row_names_intermediate=[entry for i, entry in enumerate(self.row_names) if i not in row_indices_to_remove]
            self.row_names=row_names_intermediate

    def _test_validity_of_custom_constraints(self):
        """
        Checks if IDs of added custom constraints are already present in RBA problem.
        Also checks if constraints to be removed are not among RBA constraints or are among added constraints.
        Throws error if any condition is true.
        """
        for i in self._blocks.custom_constraints.ids:
            if i in self.row_names:
                print('ID of added custom constraint {} already among constraint IDs. Please give unique name'.format(i))
                raise Exception('Invalid custom constraint.')
        for i in self._blocks.custom_constraints.constraints_to_remove:
            if i not in self.row_names:
                print('ID of constraint to remove {} not among original constraint IDs.'.format(i))
                raise Exception('Invalid custom constraint.')
            if i in self._blocks.custom_constraints.ids:
                print('ID of constraint to remove {} among custom constraints to add.'.format(i))
                raise Exception('Invalid custom constraint.')

def _build_metabolite_constraints(stoichiomerty_matrix_metabolism,
                                 enzymes,
                                 processes,
                                 targets_metabolite_components,
                                 compartments,
                                 growth_rate):
    """
    Build constraints, representing metabolite mass balance under growth and macromolecule degradation.
    Rows represent individual metabolites and columns their production/consumption 
    (metabolic reaction-stoichiometries from stoichiomerty_matrix_metabolism, flux between metabolic-network enzyme- and process-machineries and also targets.)
    """
    metabolite_constraits_growth_rate_dependent = hstack([stoichiomerty_matrix_metabolism,
                                                          growth_rate * enzymes.machinery.composition,
                                                          growth_rate * processes.machinery.composition,
                                                          targets_metabolite_components,
                                                          compartments.compute_compartment_composition_metabolite_rows(growth_rate)
                                                          ])
    
    if any(i!=0.0 for i in list(list(enzymes.machinery.complex_degradation_rate)+list(processes.machinery.complex_degradation_rate))):

        machinery_decay_rates=coo_matrix(numpy.array([list([0.0 for i in range(stoichiomerty_matrix_metabolism.shape[1])]+
                                                    list(enzymes.machinery.complex_degradation_rate)+ #enzyme.machinery.complex_hl
                                                    list(processes.machinery.complex_degradation_rate)+
                                                    [0.0]*targets_metabolite_components.shape[1]+
                                                    [0.0]*len(compartments.compartment_ids))]*stoichiomerty_matrix_metabolism.shape[0]))

        metabolite_constraits_machinery_decay_dependent = hstack([stoichiomerty_matrix_metabolism,
                                                                enzymes.machinery.composition,
                                                                processes.machinery.composition,
                                                                targets_metabolite_components,
                                                                coo_matrix((stoichiomerty_matrix_metabolism.shape[0], len(compartments.compartment_ids)))]).multiply(machinery_decay_rates)

        metabolite_constraits_degradation_machinery_decay_dependent = hstack([stoichiomerty_matrix_metabolism,
                                                                              enzymes.machinery.degradation_composition,
                                                                              processes.machinery.degradation_composition,
                                                                              targets_metabolite_components,
                                                                              coo_matrix((stoichiomerty_matrix_metabolism.shape[0], len(compartments.compartment_ids)))]).multiply(machinery_decay_rates)

        return(coo_matrix(metabolite_constraits_growth_rate_dependent.toarray()+metabolite_constraits_machinery_decay_dependent.toarray()+metabolite_constraits_degradation_machinery_decay_dependent.toarray()))

    else:

        return(metabolite_constraits_growth_rate_dependent)

def _build_process_capacity_constraints(empty_process_reaction_block,
                                       enzymes,
                                       processes,
                                       targets_machinery_cost,
                                       process_capacities,
                                       compartments,
                                       growth_rate):
    """
    Build constraints, representing process demand and capacity under growth and macromolecule degradation.
    Rows represent individual processes and columns their demand of various machineries and targets. 
    """

    process_capacity_constraints_growth_rate_dependent = hstack([empty_process_reaction_block,
                                                                            growth_rate * enzymes.machinery.processing_cost,
                                                                            growth_rate * processes.machinery.processing_cost - diags(process_capacities),
                                                                            targets_machinery_cost,
                                                                            compartments.compute_compartment_composition_process_rows(growth_rate)])

    if any(i!=0.0 for i in list(list(enzymes.machinery.complex_degradation_rate)+list(processes.machinery.complex_degradation_rate))):
        machinery_decay_rates=coo_matrix(numpy.array([list([0.0 for i in range(empty_process_reaction_block.shape[1])]+
                                                    list(enzymes.machinery.complex_degradation_rate)+
                                                    list(processes.machinery.complex_degradation_rate)+
                                                    [0.0]*targets_machinery_cost.shape[1]+
                                                    [0.0]*len(compartments.compartment_ids))]*empty_process_reaction_block.shape[0]))        
        process_capacity_constraints_machinery_decay_dependent = hstack([empty_process_reaction_block,
                                                                         enzymes.machinery.processing_cost,
                                                                         processes.machinery.processing_cost,
                                                                         targets_machinery_cost,
                                                                         coo_matrix((empty_process_reaction_block.shape[0], len(compartments.compartment_ids)))]).multiply(machinery_decay_rates)
        process_capacity_constraints_degradation_machinery_decay_dependent = hstack([empty_process_reaction_block,
                                                                                  enzymes.machinery.degradation_processing_cost,
                                                                                  processes.machinery.degradation_processing_cost,
                                                                                  targets_machinery_cost,
                                                                                  coo_matrix((empty_process_reaction_block.shape[0], len(compartments.compartment_ids)))]).multiply(machinery_decay_rates)
        
        return(coo_matrix(process_capacity_constraints_growth_rate_dependent.toarray()+process_capacity_constraints_machinery_decay_dependent.toarray()+process_capacity_constraints_degradation_machinery_decay_dependent.toarray()))
    
    else:
    
        return(process_capacity_constraints_growth_rate_dependent)

