"""Module with parameter container class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
from rba.core.functions import build_function, build_aggregate, mu_equals_mu_function, zero_function , one_function , default_ub

class Parameters(object):
    """
    Class storing RBA parameters (including aggregates).

    Attributes
    ----------
    parameters : dict
        Mapping from parameter id with object used to compute it.

    """

    def __init__(self, functions, aggregates):
        """
        Constructor.

        Parameters
        ----------
        functions : rba.xml.ListOfFunctions
            structure containing function information.
        aggregates: rba.xml.ListOfAggregates
            structure containing aggregate information.

        """
        self.parameters = {}
        self._growth_rate_fn = []
        self._medium_fn = []
        self._growth_rate_agg = []
        self._medium_agg = []
        
        # Build function parameters:
        for fn in functions:
            params = {p.id: p.value for p in fn.parameters}
            new_fn = build_function(fn.type, params, fn.variable)
            self.parameters[fn.id] = new_fn
            if new_fn.is_growth_rate_dependent() or new_fn.is_growth_rate_and_medium_dependent():
                self._growth_rate_fn.append(new_fn)
            if new_fn.is_medium_dependent() or new_fn.is_growth_rate_and_medium_dependent():
                self._medium_fn.append(new_fn)
        # add default functions to parameter set for internal use only.
        # For practical purposes in instances of Constraint class.
        self.parameters["default_function_CONSTANT_ZERO"]=zero_function
        self.parameters["default_function_CONSTANT_ONE"]=one_function
        self.parameters["default_function_GROWTHRATE"]=mu_equals_mu_function
        self.parameters["default_UB"]=default_ub
        self._growth_rate_fn.append(mu_equals_mu_function)

        # Build aggregate parameters:
        self._test_aggregate_validity(aggregates,functions)
        aggregates_accounted=0
        while aggregates_accounted < len(aggregates):
            for agg in [i for i in aggregates if i.id not in self.parameters]:
                if len([agg_ref.aggregate for agg_ref in agg.aggregate_references if agg_ref.aggregate not in self.parameters])==0:
                    new_agg = build_aggregate(agg=agg, known_functions_and_aggregates=self.parameters)
                    self.parameters[agg.id] = new_agg
                    if new_agg.is_growth_rate_dependent():
                        self._growth_rate_agg.append(new_agg)
                    if new_agg.is_medium_dependent():
                        self._medium_agg.append(new_agg)
                    aggregates_accounted+=1


    def __getitem__(self, parameter_id):
        """
        Get function or aggregate matching given id.

        Parameters
        ----------
        fn_id : str
            id of parameter to retrieve.

        Returns
        -------
        Function object.

        """
        return self.parameters[parameter_id]

    def update_growth_rate(self, growth_rate, medium):
        """
        Compute parameters for given growth rate.

        Parameters
        ----------
        growth_rate : float
            current growth rate (default 0)
        medium : dict
            Mapping of metabolite prefixes with their concentration.
        """
        # update functions first, then aggregates !!!
        for fn in self._growth_rate_fn:
            function_input=[growth_rate]
            if fn in self._medium_fn:
                variables=fn.variable.split(",")
                for variable in variables:
                    if variable == "growth_rate":
                        continue
                    elif variable in medium:
                        function_input.append(medium[variable])
                    elif variable.rsplit('_', 1)[0] in medium:
                        function_input.append(medium[variable.rsplit('_', 1)[0]])
            fn.update(*tuple(function_input))
        for agg in self._growth_rate_agg:
            agg.update()

    def update_medium(self,medium, growth_rate=0.0):
        """
        Compute parameters for given medium concentrations and growth rate.

        Parameters
        ----------
        medium : dict
            Mapping of metabolite prefixes with their concentration.
        growth_rate : float
            current growth rate (default 0)
        """
        # update functions first, then aggregates !!!
        for fn in self._medium_fn:
            # /!\ we identify metabolites by their prefix !!!
            variables=fn.variable.split(",")
            function_input=[]
            for variable in variables:
                if variable == "growth_rate":
                    function_input.append(growth_rate)
                elif variable in medium:
                    function_input.append(medium[variable])
                elif variable.rsplit('_', 1)[0] in medium:
                    function_input.append(medium[variable.rsplit('_', 1)[0]])
            fn.update(*tuple(function_input))
        for agg in self._medium_agg:
            agg.update()

    def _test_aggregate_validity(self,aggregates,functions):
        """
        Tests if all referenced operands in all aggregates are among defined parameters.
        If not an error is raised.
        """
        function_ids=[i.id for i in functions]
        aggregate_ids=[i.id for i in aggregates]
        for agg in aggregates:
            if agg.id in function_ids:
                    print('Aggregate ID {} already among function IDs '.format(agg.id))
                    raise Exception('Invalid aggregate.')
            for function_ref in agg.function_references:
                if function_ref.function not in function_ids:
                    print('Unknown function reference {} in aggregate {} '.format(function_ref.function,agg.id))
                    raise Exception('Invalid aggregate.')
            for agg_ref in agg.aggregate_references:
                if agg_ref.aggregate not in aggregate_ids:
                    print('Unknown aggregate reference {} in aggregate {} '.format(agg_ref.aggregate,agg.id))
                    raise Exception('Invalid aggregate.')
