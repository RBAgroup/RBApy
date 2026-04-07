"""Module containing valid RBA functions."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import numpy


class BaseFunction(object):
    """Mother of all functions."""

    def __init__(self, variable):
        """Build default object."""
        self.variable = variable

    def is_growth_rate_dependent(self):
        """Return whether function depends on growth rate."""
        return self.variable == 'growth_rate'

    def is_medium_dependent(self):
        """Return whether function depends on medium variable."""
        return self.variable and ('growth_rate' not in self.variable)

    def is_growth_rate_and_medium_dependent(self):
        """Return whether function depends on growth rate."""
        return self.variable and (self.variable != 'growth_rate') and ('growth_rate' in self.variable)

class Keff_psi2(BaseFunction):
    """
    Apparent catalytic rate of photosystem 2 for Arabidopsis thaliana.
    """

    name = 'psi2'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            OBSOLETE 
        variable : str
            Function variable(s).

        """
        super(Keff_psi2, self).__init__(variable)
        self.value = 0

    def update(self, CO2,O2,Temperature,hnu):
        """Evaluate function."""
        self.value = rubisco_psii_update(T_celcius=Temperature,CO2=CO2,O2=O2,hnu=hnu)["keff_psiII"]

    def is_medium_dependent(self):
        """Return whether function depends on medium variable."""
        return True


class Keff_Rubisco_Oxy(BaseFunction):
    """
    Apparent catalytic rate of Rubisco-associated photorespiration in Arabidopsis thaliana.
    """

    name = 'rubiscoOxygenase'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            OBSOLETE 
        variable : str
            Function variable(s).

        """
        super(Keff_Rubisco_Oxy, self).__init__(variable)
        self.value = 0

    def update(self, CO2,O2,Temperature,hnu):
        """Evaluate function."""
        self.value = rubisco_psii_update(T_celcius=Temperature,CO2=CO2,O2=O2,hnu=hnu)["keff_rubisco_oxygenase"]
        
    def is_medium_dependent(self):
        """Return whether function depends on medium variable."""
        return True


class Keff_Rubisco_Carboxy(BaseFunction):
    """
    Apparent catalytic rate of Rubisco-associated carbon fixation in Arabidopsis thaliana.
    """

    name = 'rubiscoCarboxygenase'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            OBSOLETE 
        variable : str
            Function variable(s).

        """
        super(Keff_Rubisco_Carboxy, self).__init__(variable)
        self.value = 0

    def update(self, CO2,O2,Temperature,hnu):
        """Evaluate function."""
        self.value = rubisco_psii_update(T_celcius=Temperature,CO2=CO2,O2=O2,hnu=hnu)["keff_rubisco_carboxylase"]

    def is_medium_dependent(self):
        """Return whether function depends on medium variable."""
        return True


class LogisticFunction(BaseFunction):
    """
    Class computing logistic functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'logistic'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: RATE.
        variable : str
            Function variable.

        """
        super(LogisticFunction, self).__init__(variable)
        self._a=parameters['A']
        self._b=parameters['B']
        self._c=parameters['C']
        self._d=parameters['D']
        self.value = self._a / (self._b + self._c)

    def update(self, x):
        """Evaluate function."""
        self.value = self._a / (self._b + (self._c * numpy.exp(-self._d * x)))


class ConstantFunction(BaseFunction):
    """
    Class computing constant functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'constant'

    def __init__(self, parameters, variable=None):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: CONSTANT.
        variable : str
            Function variable.

        """
        super(ConstantFunction, self).__init__(None)
        self.value = parameters['CONSTANT']

    def update(self, x):
        """Evaluate function."""
        pass


class ExponentialFunction(BaseFunction):
    """
    Class computing exponential functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'exponential'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: RATE.
        variable : str
            Function variable.

        """
        super(ExponentialFunction, self).__init__(variable)
        self._multiplier = parameters.get('MULTIPLIER',1.0)
        self._rate = parameters['RATE']
        self._constant = parameters.get('CONSTANT', 0.0)
        self.value = self._multiplier

    def update(self, x):
        """Evaluate function."""
        self.value = self._multiplier*numpy.exp(self._constant + self._rate * x)


class IndicatorFunction(BaseFunction):
    """
    Class computing indicator functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : bool
        current function value.

    """

    name = 'indicator'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: X_MIN, X_MAX.
        variable : str
            Function variable.

        """
        super(IndicatorFunction, self).__init__(variable)
        self.value = 0
        self._x_min = parameters['X_MIN']
        self._x_max = parameters['X_MAX']

    def update(self, x):
        """Evaluate function."""
        self.value = float((x > self._x_min) and (x < self._x_max))


class LinearFunction(BaseFunction):
    """
    Class computing linear functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'linear'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: X_MIN, X_MAX,
            LINEAR_COEF, LINEAR_CONSTANT, Y_MIN, Y_MAX.
        variable : str
            Function variable.

        """
        super(LinearFunction, self).__init__(variable)
        self._x_min = parameters['X_MIN']
        self._x_max = parameters['X_MAX']
        self._coef = parameters['LINEAR_COEF']
        self._constant = parameters['LINEAR_CONSTANT']
        self._y_min = parameters['Y_MIN']
        self._y_max = parameters['Y_MAX']
        self.value = float(min(max(self._constant, self._y_min), self._y_max))

    def update(self, x):
        """Evaluate function."""
        x_eval = min(max(x, self._x_min), self._x_max)
        y = self._coef * x_eval + self._constant
        self.value = float(min(max(y, self._y_min), self._y_max))


class QuadraticFunction(BaseFunction):
    """
    Class computing quadratic functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'quadratic'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: X_MIN, X_MAX,
            QUADRATIC_TERM_COEF, LINEAR_TERM_COEF, CONSTANT, Y_MIN, Y_MAX.
        variable : str
            Function variable.

        """
        super(LinearFunction, self).__init__(variable)
        self._x_min = parameters['X_MIN']
        self._x_max = parameters['X_MAX']
        self._coef_quadratic_term = parameters['QUADRATIC_TERM_COEF']
        self._coef_linear_term = parameters['LINEAR_TERM_COEF']
        self._constant = parameters['CONSTANT']
        self._y_min = parameters['Y_MIN']
        self._y_max = parameters['Y_MAX']
        self.value = float(min(max(self._constant, self._y_min), self._y_max))

    def update(self, x):
        """Evaluate function."""
        x_eval = min(max(x, self._x_min), self._x_max)
        y = self._coef_quadratic_term * x_eval**2 + self._coef_linear_term * x_eval + self._constant
        self.value = min(max(y, self._y_min), self._y_max)


class MichaelisMentenFunction(BaseFunction):
    """
    Class computing michaelis menten functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'michaelisMenten'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: kmax, Km.
            Optionally, it may contain Y_MIN and HILL_COEFFICIENT (1 by default).
        variable : str
            Function variable.

        """
        super(MichaelisMentenFunction, self).__init__(variable)
        self._kmax = float(parameters['kmax'])
        self._km = float(parameters['Km'])
        self._y_min = parameters.get('Y_MIN', None)
        self._n = parameters.get('HILL_COEFFICIENT', 1.0)
        self.value = 0

    def update(self, x):
        """Evaluate function at given point."""
        if x !=0:
            y = self._kmax * ((x**self._n) / ((x**self._n) + (self._km**self._n)))
        else:
            y=0
        self.value = max(y, self._y_min) if self._y_min else y


class CompetitiveInhibitionFunction(BaseFunction):
    """
    Class computing michaelis menten functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'competitiveInhibition'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: kmax, Km, Ki.
            Optionally, it may contain Y_MIN.
        variable : str
            Function variable.

        """
        super(CompetitiveInhibitionFunction, self).__init__(variable)
        self._kmax = float(parameters['kmax'])
        self._Km = float(parameters['Km'])
        self._Ki = float(parameters['Ki'])
        self._y_min = parameters.get('Y_MIN', None)
        self.value = 0

    def update(self, x_1, x_2):
        """Evaluate function at given point."""
        y = self._kmax * x_1 / (x_1 + self._Km * (1 + x_2 / self._Ki))
        self.value = max(y, self._y_min) if self._y_min else y


class InverseFunction(BaseFunction):
    """
    Class computing inverse functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'inverse'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: CONSTANT.
        variable : str
            Function variable.

        """
        super(InverseFunction, self).__init__(variable)
        self._xmax = parameters['CONSTANT']
        self.value = 0

    def update(self, x):
        """Evaluate function."""
        try:
            self.value = self._xmax / x
        except KeyError:
            print('variable is 0, impossible to do inversion')


class ArrheniusFunction(BaseFunction):
    """
    Class computing arrhenius functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'arrhenius'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: Y_MIN, Y_MAX,
            PRE_EXPONENTIAL_FACTOR, ACTIVATION_ENERGY and GAS_CONSTANT.
        variable : str
            Function variable.

        """
        super(ArrheniusFunction, self).__init__(variable)
        self.value = 0
        self._y_min = parameters['Y_MIN']
        self._y_max = parameters['Y_MAX']
        self._activation_energy = parameters['ACTIVATION_ENERGY']
        self._A = parameters['PRE_EXPONENTIAL_FACTOR']
        self._R = parameters['GAS_CONSTANT']

    def update(self, x):
        """Evaluate function at given point."""
        try:
            y=self._A*numpy.exp(-(self._activation_energy/(self._R * (x+273))))
            self.value = min(max(y, self._y_min), self._y_max)
        except KeyError:
            print('variable is 0, impossible to do inversion')


class CenteredArrheniusFunction(BaseFunction):
    """
    Class computing centered arrhenius functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    variable : str
        variable used by function.
    value : float
        current function value.

    """

    name = 'centeredArrhenius'

    def __init__(self, parameters, variable):
        """
        Constructor.

        Parameters
        ----------
        parameters : dict
            Dict that must contain the following keys: Y_MIN, Y_MAX,
            REFERENCE_TEMPERATURE, PRE_EXPONENTIAL_FACTOR, ACTIVATION_ENERGY and GAS_CONSTANT.
        variable : str
            Function variable.

        """
        super(CenteredArrheniusFunction, self).__init__(variable)
        self.value = parameters['PRE_EXPONENTIAL_FACTOR']
        self._y_min = parameters['Y_MIN']
        self._y_max = parameters['Y_MAX']
        self._activation_energy = parameters['ACTIVATION_ENERGY']
        self._A = parameters['PRE_EXPONENTIAL_FACTOR']
        self._R = parameters['GAS_CONSTANT']
        self._Tref = parameters['REFERENCE_TEMPERATURE']

    def update(self, x):
        """Evaluate function at given point."""
        try:
            y=self._A*numpy.exp(-(((1/(x+273))-(1/(self._Tref+273)))*(self._activation_energy/self._R)))
            self.value = min(max(y, self._y_min), self._y_max)
        except KeyError:
            print('variable is 0, impossible to do inversion')


class BaseAggregate(object):
    """Mother of all aggregates."""

    def __init__(self, operand_handles, operand_exponents):
        """Build default object."""
        self._operands = operand_handles
        self._operand_exponents = operand_exponents

    def is_growth_rate_and_medium_dependent(self):
        """Return whether aggregate depends on growth rate."""
        return any([op.is_growth_rate_and_medium_dependent() for op in self._operands])

    def is_growth_rate_dependent(self):
        """Return whether aggregate depends on growth rate."""
        return any([op.is_growth_rate_dependent() for op in self._operands]+[op.is_growth_rate_and_medium_dependent() for op in self._operands])

    def is_medium_dependent(self):
        """Return whether aggregate depends on medium variable."""
        return any([op.is_medium_dependent() for op in self._operands]+[op.is_growth_rate_and_medium_dependent() for op in self._operands])


class MultiplicationAggregate(BaseAggregate):
    """
    Class computing multiplication of operands.

    Attributes
    ----------
    name : str
        identifier of this class.
    value : int
        current function value.

    """

    name = 'multiplication'

    def __init__(self, operand_handles, operand_exponents):
        """
        Constructor.

        Parameters
        ----------
        operand_handles : list of operand handles
            Operands that have to be multiplied.
        operand_exponents : list of operand exponents
            Exponents of operands that have to be multiplied.

        """
        super(MultiplicationAggregate, self).__init__(operand_handles, operand_exponents)
        self.value = 1
        self.update()


    def update(self):
        """Evaluate function."""
        y = 1
        for i in range(len(self._operands)):
            if self._operands[i].value !=0:
                y *= (self._operands[i].value)**self._operand_exponents[i]
            else:
                y *= 0
        self.value = y


class AdditionAggregate(BaseAggregate):
    """
    Class computing addition of functions.

    Attributes
    ----------
    name : str
        identifier of this class.
    value : int
        current function value.

    """

    name = 'addition'

    def __init__(self, operand_handles, operand_exponents):
        """
        Constructor.

        Parameters
        ----------
        operand_handles : list of operand handles
            Operands that have to be added up.
        operand_exponents : list of operand exponents
            Exponents of operands that have to be multiplied.

        """
        super(AdditionAggregate, self).__init__(operand_handles, operand_exponents)
        self.value = 0
        self.update()


    def update(self):
        """Evaluate function."""
        y = 0
        for i in range(len(self._operands)):
            if self._operands[i].value !=0:
                y += (self._operands[i].value)**self._operand_exponents[i]
        self.value = y


# list of accepted function names and classes implementing them
SIMPLE = [ConstantFunction, LinearFunction, QuadraticFunction, IndicatorFunction,
          ExponentialFunction, MichaelisMentenFunction, InverseFunction,
          CompetitiveInhibitionFunction,ArrheniusFunction,CenteredArrheniusFunction,
          Keff_psi2, Keff_Rubisco_Oxy, Keff_Rubisco_Carboxy, LogisticFunction]
AGGREGATE = [MultiplicationAggregate, AdditionAggregate]
VALID_FNS = {c.name: c for c in SIMPLE}
VALID_AGGS = {c.name: c for c in AGGREGATE}


def build_function(type_, params, variable):
    """
    Create object function matching type and parameters given.

    Parameters
    ----------
    type_ : str
        'name' attribute of one of the classes of this module.
    params : dict
        Dict mapping parameter names with their values.
    variable : str
        Function variable.

    Returns
    -------
    Function object matching type and parameters provided.

    """
    try:
        # retrieve class implementing function
        fn_class = VALID_FNS[type_]
    except KeyError:
        print('Unknown function type: ' + type_ + '. Valid types are: '
              + ', '.join(VALID_FNS.keys()))
        raise UserWarning('Invalid function.')
    try:
        return fn_class(params, variable)
    except KeyError as error:
        print('Missing parameter: ' + str(error))
        raise UserWarning('Invalid function.')


def build_aggregate(agg, known_functions_and_aggregates):
    """
    Create aggregate from xml_structure.

    Parameters
    ----------
    agg : XML node
        Structure describing aggregate.
    known_functions_and_aggregates : rba.core.parameters.Parameters
        Known parameters.

    Returns
    -------
    Aggregate object matching xml structure.

    """
    try:
        # retrieve class implementing aggregate
        agg_class = VALID_AGGS[agg.type]
    except KeyError:
        print('Unknown aggregate type: ' + agg.type + '. Valid types are: '
              + ', '.join(VALID_AGGS.keys()))
        raise UserWarning('Invalid aggregate.')
    
    try:
        # retrieve functions used in aggregate
        fn_handles = [known_functions_and_aggregates[ref.function] for ref in agg.function_references]
        fn_exponents= [ref.exponent for ref in agg.function_references]
    except KeyError as error:
        print('Unknown function reference: ' + error.args[0])
        raise UserWarning('Invalid aggregate.')
    
    try:
        # retrieve aggregates used in aggregate
        agg_handles = [known_functions_and_aggregates[ref.aggregate] for ref in agg.aggregate_references]
        agg_exponents= [ref.exponent for ref in agg.aggregate_references]
    except KeyError as error:
        print('Unknown aggregate reference: ' + error.args[0])
        raise UserWarning('Invalid aggregate.')
    
    # initiate aggregate with referenced functions and aggregates as operands
    return agg_class(list(fn_handles+agg_handles),list(fn_exponents+agg_exponents))

# define functions to be imported in other modules
default_lb = ConstantFunction({'CONSTANT': -1e3})
default_ub = ConstantFunction({'CONSTANT': 1e5})
zero_function=ConstantFunction({'CONSTANT': 0.0})
one_function=ConstantFunction({'CONSTANT': 1.0})
mu_equals_mu_function=LinearFunction(parameters={'LINEAR_COEF': 1.0, 'LINEAR_CONSTANT': 0.0,'X_MIN': 0, 'X_MAX': float('Inf'),'Y_MIN': 0.0, 'Y_MAX': float('Inf')},variable='growth_rate')

def rubisco_psii_update(T_celcius,CO2,O2,hnu):
    # Values for kcat, Kc, Km at 25°C were extracted 
    # from  walker et al. Plant cell env, 2014.
    T_kelvin = T_celcius + 273.15
    R = 8.31 # J/K/mol
    T_sat_rubisco = 310 # kelvin
    aR = 1/(1+numpy.exp(0.3*(T_kelvin-T_sat_rubisco))) # activation rubisco
    # From Biochemical Model of C 3 Photosynthesis
    # von caemmer et al. 
    E_Ko = 35.94e3 #J/mol
    E_Kc = 59.36e3 #J/mol
    E_Vcmax = 58.52e3 # J/mol
    # from table 5 walker et al. Plant cell env, 2014.
    Kc_25 = 312e-3 # mbar
    Ko_25 =  181 # mbar
    ratio_Vomax_Vc_max_25 = 0.24
    # from table 3 walker et al. Plant cell env, 2014.
    kcat_c_25 = 3.3*3600*24*1.5
    kcat_o_25 = kcat_c_25*ratio_Vomax_Vc_max_25
    # calcul kcat_c, kcat_o, Kc, Ko à la temperature d'interet
    kcat_c =aR*kcat_c_25*numpy.exp((T_celcius-25)*E_Vcmax/(298*R*T_kelvin)) # 4.1 1/s
    kcat_o = aR*kcat_o_25*numpy.exp((T_celcius-25)*E_Vcmax/(298*R*T_kelvin)) # 1 1/s
    Kc = Kc_25/1.01325*numpy.exp((T_celcius-25)*E_Kc/(298*R*T_kelvin)) # 9.8µM converti en mM
    Ko = Ko_25/1.01325*numpy.exp((T_celcius-25)*E_Ko/(298*R*T_kelvin)) # 470µM
    # ratio Jmax/Vcmax
    # calcule from Table  4 of walker et al. Plant cell env, 2014.
    # normalement entre 1.5 et 2 à 25°C. 
    # 2.44 ds walker et al .
    ratio_Jmax25_Vcmax25 = 2.4441
    kcat_psiII_25 = kcat_c_25*ratio_Jmax25_Vcmax25*(hnu/(1+hnu)) 
    E_Jmax = 37e3 # J/mol
    H = 220e3 #J/mol
    S = 710 # J/K/mol
    keff_psiII = kcat_psiII_25*numpy.exp((T_kelvin-298)*E_Jmax/(298*R*T_kelvin))*(1+numpy.exp((298*S-H)/298/R))/(1+numpy.exp((S*T_kelvin-H)/R/T_kelvin))
    keff_rubisco_carboxylase = kcat_c/(1+ Kc/CO2 + Kc*O2/Ko/CO2)
    keff_rubisco_oxygenase = kcat_o/(1 + Ko/O2 +  Ko*CO2/Kc/O2)
    rubisco_psii_keff = {'keff_psiII': float(keff_psiII), 'keff_rubisco_carboxylase': float(keff_rubisco_carboxylase),
     'keff_rubisco_oxygenase':  float(keff_rubisco_oxygenase)}
    return(rubisco_psii_keff)