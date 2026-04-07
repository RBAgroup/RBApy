"""Module defining parameter-specific classes for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml._rbaml_version import __rbaml_version__ as rbaml_version
from rba.xml.common import get_unique_child, ListOf, xml_input_tag_error

__all__ = ['RbaParameters', 'Parameter', 'ListOfParameters', 'Function',
           'ListOfFunctions', 'FunctionReference', 'ListOfFunctionReferences',
           'Aggregate', 'ListOfAggregates', 'AggregateReference', 'ListOfAggregateReferences']


class RbaParameters(object):
    """
    Parameter-related structures.

    Attributes
    ----------
    target_densities : ListOfTargetDensities
        List of density constraints.
    functions : ListOfFunctions
        List of user-defined functions.
    aggregates : ListOfAggregated
        List of user-defined aggregates (composition of functions).

    """
    tag='RBAParameters'
    def __init__(self):
        """Constructor."""
        self.functions = ListOfFunctions()
        self.aggregates = ListOfAggregates()

    def write(self, output_stream,meta_data={}):
        """
        Write information as an XML structure.

        Parameters
        ----------
        output_stream : file or buffer
            Location where XML structure should be written.
        """
        root = etree.Element(self.tag)
        root.extend([self.functions.to_xml_node(),
                     self.aggregates.to_xml_node()])
        tree=etree.ElementTree(root)
        root.addprevious(etree.Comment(" Created by RBApy version {} with RBAML version {}".format(meta_data.get("RBApy_version_model_generation",None),rbaml_version)))
        root.addprevious(etree.Comment("Model: {}, Organism/Tissue: {}, Taxon-ID: {}".format(meta_data.get("Model",""),
                                                                                             meta_data.get("Organism/Tissue",""),
                                                                                             meta_data.get("Taxon_ID",""))))
        root.addprevious(etree.Comment("Authored by: {}".format(meta_data.get("Author(s)",""))))
        root.addprevious(etree.Comment("Copyright of: {}".format(meta_data.get("Copyright_holder",""))))
        root.addprevious(etree.Comment("License: {}".format(meta_data.get("License",""))))
        root.addprevious(etree.ProcessingInstruction("rbaml", 'version="{}"'.format(rbaml_version)))
        tree.write(output_stream, xml_declaration=True, encoding="UTF-8", pretty_print=True)

    @classmethod
    def from_file(cls, input_stream):
        """
        Build object from XML structure.

        Parameters
        ----------
        input_stream : file or buffer
            Location containing XML structure.

        """
        node = etree.ElementTree(file=input_stream).getroot()
        # Verify that root tag in file is the one specified by class RBAParameters
        if node.tag != cls.tag:
            raise xml_input_tag_error(file=input_stream.name,file_tag=node.tag,requested_tag=cls.tag)

        result = cls()
        n = get_unique_child(node, ListOfFunctions.tag)
        result.functions = ListOfFunctions.from_xml_node(n)
        n = get_unique_child(node, ListOfAggregates.tag)
        result.aggregates = ListOfAggregates.from_xml_node(n)
        return result


class Parameter(object):
    """
    Parameter represented with an id, value couple.

    Attributes
    ----------
    id : str
        Identifier of parameter.
    value : float
        Value of parameter.

    """

    tag = 'parameter'

    def __init__(self, id_, value):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier of parameter.
        value : float
            Value of parameter.

        """
        self.id = id_
        self.value = value

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('value', str(self.value))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('id'), float(node.get('value')))


class ListOfParameters(ListOf):
    """List of Parameter elements."""

    tag = 'listOfParameters'
    list_element = Parameter


class Function(object):
    """
    Function defined by a type and parameters.

    Attributes
    ----------
    id : str or None
        Identifier of function (if applicable).
    type : str
        Type of function.
    parameters : ListOfParameters
        List of parameters used by function.
    variable : str
        Name of variable (if applicable).

    """

    tag = 'function'

    def __init__(self, id_, type_, parameters=None, variable=None):
        """
        Constructor.

        Parameters
        ----------
        id_ : str or None
            identifier of function (if applicable).
        type_: str
            type of function.
        parameters : dict, optional
            dict containing parameters used by function.
        variable : str, optional
            name of variable (set to growth_rate if empty).

        """
        self.id = id_ if id_ is not None else ''
        self.type = type_
        if variable:
            self.variable = variable
        else:
            self.variable = 'growth_rate'
        self.set_parameters(parameters)

    def set_parameters(self, parameters):
        """Create parameter list from dictionary."""
        self.parameters = ListOfParameters()
        if parameters:
            for key, value in parameters.items():
                self.parameters.append(Parameter(key, value))

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        if self.variable is not None:
            result.set('variable', self.variable)
        if not self.parameters.is_empty():
            result.extend([self.parameters.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'), node.get('type'),
                     {}, node.get('variable'))
        n = get_unique_child(node, 'listOfParameters', False)
        if n is not None:
            result.parameters = ListOfParameters.from_xml_node(n)
        return result


class ListOfFunctions(ListOf):
    """List of Function elements."""

    tag = 'listOfFunctions'
    list_element = Function


class FunctionReference(object):
    """
    Reference to a Function.

    Attributes
    ----------
    function : str
        Function identifier.
    exponent : float
        Exponent, related to function in reference

    """

    tag = 'functionReference'

    def __init__(self, function,exponent="1.0"):
        """
        Constructor.

        Parameters
        ----------
        function : str
            Function identifier.
        exponent : float
            Exponent, related to function in reference

        """
        self.function = function
        self.exponent = float(exponent)

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('function', self.function)
        result.set('exponent', str(self.exponent))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        if node.get('exponent'):
            return cls(node.get('function'),node.get('exponent'))
        else:
            return cls(node.get('function'),"1.0")


class AggregateReference(object):
    """
    Reference to a Aggregate.

    Attributes
    ----------
    aggregate : str
        Aggregate identifier.
    exponent : float
        Exponent, related to aggregate in reference

    """

    tag = 'aggregateReference'

    def __init__(self, aggregate,exponent="1.0"):
        """
        Constructor.

        Parameters
        ----------
        aggregate : str
            Aggregate identifier.
        exponent : float
            Exponent, related to function in reference

        """
        self.aggregate = aggregate
        self.exponent = float(exponent)

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('aggregate', self.aggregate)
        result.set('exponent', str(self.exponent))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        if node.get('exponent'):
            return cls(node.get('aggregate'),node.get('exponent'))
        else:
            return cls(node.get('aggregate'),"1.0")


class ListOfFunctionReferences(ListOf):
    """List of FunctionReference elements."""

    tag = 'listOfFunctionReferences'
    list_element = FunctionReference


class ListOfAggregateReferences(ListOf):
    """List of AggregateReference elements."""

    tag = 'listOfAggregateReferences'
    list_element = AggregateReference


class Aggregate(object):
    """
    Aggregate (composition of Functions).

    Attributes
    ----------
    id : str
        Identifier.
    type : str
        Type of aggregation (e.g. 'multiplication').
    function_references : ListOfFunctionReferences
        References to functions aggregated.

    """

    tag = 'aggregate'

    def __init__(self, id_, type_):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        type_ : str
            Type of aggregation ('multiplication' or 'addition').

        """
        self.id = id_ if id_ is not None else ''
        self.type = type_
        self.function_references = ListOfFunctionReferences()
        self.aggregate_references = ListOfAggregateReferences()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('type', self.type)
        opt = [self.function_references, self.aggregate_references]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls(node.get('id'), node.get('type'))
        n_function_refs = get_unique_child(node, ListOfFunctionReferences.tag, strict=False)
        if n_function_refs is not None:
            result.function_references = ListOfFunctionReferences.from_xml_node(n_function_refs)
        n_aggregate_refs = get_unique_child(node, ListOfAggregateReferences.tag, strict=False)
        if n_aggregate_refs is not None:
            result.aggregate_references = ListOfAggregateReferences.from_xml_node(n_aggregate_refs)
        if n_function_refs is None and n_aggregate_refs is None:
            print('No functions or aggregates referenced in aggregate {} (at least one of the two required)'.format(node.get('id')))
            raise Exception('Invalid aggregate.')
        if len(result.function_references)==0 and len(result.aggregate_references)==0:
            print('Empty lists of function- and aggregate references in aggregate {} (At least one of the two must be non-empty)'.format(node.get('id')))
            raise Exception('Invalid aggregate.')
        return result


class ListOfAggregates(ListOf):
    """List of Aggregate elements."""

    tag = 'listOfAggregates'
    list_element = Aggregate
