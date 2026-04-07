"""Module defining CustomConstraint-related classes."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy
from scipy.sparse import coo_matrix
# local imports
from rba.core.parameter_vector import ParameterVector


class CustomConstraints(object):
    """
    Class computing substructures, representing custom linear constraints imposed on model.

    Constraint ID: sum(a_1*v_1 , a_2*v_2 , a_3*v_3 , ...) =/<=/>= c - sum(Param_1(µ)*x_1 , Param_2(µ)*x_2 , ...)

     Attributes
    ----------
    ids :  list of str
        IDs of custom constraints.
        
    signs : list of str
        Signs of custom constraints (E: = , L: <= and G: >=).

    variable_dependent_lhs_variables : list of lists of str
        IDs of decision variables in RBA-problem in linear sum on LHS (v_...)
        
    variable_dependent_lhs_constant_coefficients : list of lists of float
        Coefficients of variable_dependent_lhs_variables in linear sum on LHS (a_...)
                
    constant_rhs :  list of float
        Constant (parameter-independent) term of righthand side (c).
        
    parameter_dependent_rhs_parameters : list of rba.core.parameter_vector.ParameterVector
        Growth-rate specific evaluated parameter-values (Param_...(µ))
        
    parameter_dependent_rhs_constant_coefficients : list of numpy.arrays
        Coefficients of parameter_dependent_rhs_parameters in linear sum on RHS (x_...)
    """

    def __init__(self, custom_constraints, parameters):
        """
        Constructor.

        Parameters
        ----------
        custom_constraints : rba.xml.RbaCustomConstraints
            Data structure with custom constraint info.
        parameters : rba.core.parameters.Parameters
            Parameter information
        """

        # initiate empty list and populate them with the entries for different constraints
        self.ids = []
        self.signs=[]
        self.constant_rhs=[]
        self.variable_dependent_lhs_variables=[]
        self.variable_dependent_lhs_constant_coefficients=[]
        self.variable_dependent_lhs_parameter_coefficients=[]
        self.variable_dependent_lhs_growth_rate_proportionality_coefficients=[]
        self.parameter_dependent_rhs_parameters=[]
        self.parameter_dependent_rhs_constant_coefficients=[]
        self.parameter_dependent_rhs_growth_rate_proportionality_coefficients=[]
        self.constraints_to_remove = []
        for constraint in custom_constraints.constraints:
            if constraint.action=="remove":
                self.constraints_to_remove.append(constraint.id)
            elif constraint.action=="add":
                self.ids.append(constraint.id)
                self.variable_dependent_lhs_variables.append([v_ref.variable for v_ref in constraint.definition.variable_references])
                self.variable_dependent_lhs_constant_coefficients.append([v_ref.constant_coefficient for v_ref in constraint.definition.variable_references])  
                self.variable_dependent_lhs_parameter_coefficients.append(ParameterVector(parameter_list=[v_ref.parameter_coefficient if v_ref.parameter_coefficient else "default_function_CONSTANT_ONE" for v_ref in constraint.definition.variable_references], known_parameters=parameters))  
                self.variable_dependent_lhs_growth_rate_proportionality_coefficients.append(ParameterVector(parameter_list=["default_function_GROWTHRATE" if v_ref.multiplied_with_growth_rate else "default_function_CONSTANT_ONE" for v_ref in constraint.definition.variable_references], known_parameters=parameters))              
                self.parameter_dependent_rhs_parameters.append(ParameterVector(parameter_list=[p_ref.parameter for p_ref in constraint.definition.parameter_references], known_parameters=parameters))
                self.parameter_dependent_rhs_constant_coefficients.append(numpy.array([p_ref.constant_coefficient for p_ref in constraint.definition.parameter_references]))            
                self.parameter_dependent_rhs_growth_rate_proportionality_coefficients.append(ParameterVector(parameter_list=["default_function_GROWTHRATE" if p_ref.multiplied_with_growth_rate else "default_function_CONSTANT_ONE" for p_ref in constraint.definition.parameter_references], known_parameters=parameters))  

                if constraint.definition.value is not None:
                    self.constant_rhs.append(float(constraint.definition.value))
                    self.signs.append('E')
                elif constraint.definition.upper_bound is not None:
                    self.constant_rhs.append(float(constraint.definition.upper_bound))
                    self.signs.append('L')
                elif constraint.definition.lower_bound is not None:
                    self.constant_rhs.append(float(constraint.definition.lower_bound))
                    self.signs.append('G')
                else:
                    raise UserWarning('Constraint ' + constraint.id
                                    + ': you must specify a value or a '
                                    'bound.')
            elif constraint.action=="ignore":
                continue
            
    def build_custom_constraints_lhs(self,col_names):
        """
        Build lefthand side of custom constraints (all custom constraints stacked).
        """
        #Initiate empty matrix with n= number of custom constraints rows, 
        #and m=number of decision variables in rba problem cols:
        custom_constraints=numpy.zeros((len(self.ids),len(col_names)))

        #Populate zero matrix with linear coefficients at the position of the respecive decision variable:
        #Iterate over all custom constraints:
        for i in range(len(self.ids)):
            row_indices=[i]*len(self.variable_dependent_lhs_variables[i])
            #find position of decision variables referenced in constraint:
            col_indices=[col_names.index(j) for j in self.variable_dependent_lhs_variables[i]]
            #replace zeros with respective coefficients:
            
            custom_constraints[row_indices,col_indices]=[const_coeff*param_coeff*mu_coeff for const_coeff, param_coeff, mu_coeff in zip(self.variable_dependent_lhs_constant_coefficients[i],self.variable_dependent_lhs_parameter_coefficients[i].compute(),self.variable_dependent_lhs_growth_rate_proportionality_coefficients[i].compute())]
        #return as sparse matrix:
        return(coo_matrix(custom_constraints))

    def build_custom_constraints_rhs(self):
        """
        Build growth rate dependent righthand side of custom constraints (all custom constraints stacked).
        """
        return(numpy.array([self.constant_rhs[i]- #constant RHS term (subtract parameter terms)
                            numpy.sum(numpy.array([const_coeff*param_coeff*mu_coeff for const_coeff, param_coeff, mu_coeff in zip(self.parameter_dependent_rhs_constant_coefficients[i], self.parameter_dependent_rhs_parameters[i].compute(),self.parameter_dependent_rhs_growth_rate_proportionality_coefficients[i].compute())]))
                            for i in range(len(self.ids))] # do for all constraints and stack
                        )
            )

