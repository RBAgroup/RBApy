"""
Module defining classes used for defining constraint RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml._rbaml_version import __rbaml_version__ as rbaml_version
from rba.xml.common import get_unique_child, ListOf, is_true, xml_input_tag_error

__all__ = ['RbaCustomConstraints', 'Constraint', 'ConstraintDefinition','VariableReference','ConstraintReference', 
           'ParameterReference', 'ListOfConstraints', 'ListOfVariableReferences','ListOfParameterReferences']

class RbaCustomConstraints(object):
    """
    Custom constraints related structures.

    Attributes
    ----------
    constraints : ListOfLinearConstraints
        List of user defined custom constraints to add.

    """
    tag='RBACustomConstraints'
    def __init__(self):
        """
        Default constructor.
        """
        self.constraints = ListOfConstraints()

    def write(self, output_stream,meta_data={}):
        """
        Write information as an XML structure.

        Parameters
        ----------
        output_stream : file or buffer
            Location where XML structure should be written.
        """
        root = etree.Element(self.tag)
        root.extend([self.constraints.to_xml_node()])
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
        Constructor from XML structure.

        Parameters
        ----------
        input_stream : file or buffer
            Location containing XML structure.
        """
        node = etree.ElementTree(file=input_stream).getroot()
        # Verify that root tag in file is the one specified by class RBACustomConstraints
        if node.tag != cls.tag:
            raise xml_input_tag_error(file=input_stream.name,file_tag=node.tag,requested_tag=cls.tag)
        result = cls()
        n = get_unique_child(node, ListOfConstraints.tag)
        result.constraints = ListOfConstraints.from_xml_node(n)
        return result


class Constraint(object):
    """
    Constraint.

    Attributes
    ----------
    id: str
        Identifier.
    action : str
        add, remove, ignore (default add)
    """
    tag = 'constraint'

    def __init__(self):
        """Build default object."""
        self.id = None
        self.action = None
        self.definition = ConstraintDefinition()
    
    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('action', self.action)
        # optional subelements
        opt = [self.definition]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls()
        result._init_from_xml_node(node)
        n = get_unique_child(node, ConstraintDefinition.tag, False)
        if n is not None:
            result.definition = ConstraintDefinition.from_xml_node(n)
        return result

    def _init_from_xml_node(self, node):
        """Match attributes with given node."""
        self.id = node.get('id')
        self.action = node.get('action')
        self.definition = ConstraintDefinition()


class ConstraintDefinition(object):
    """
    ConstraintDefinition.

    Attributes
    ----------
    value : float or None
        exact target value (None if no exact value to match).
    lower_bound : float or None
        lower bound on target value (None if no lower bound).
    upper_bound : float or None
        upper bound on target value (None if no upper bound).
    variable_references : ListOfVariableReferences
        Variable references in constraint
    parameter_references : ListOfParameterReferences
        Parameter references in constraint
    """

    tag = 'constraintDefinition'

    def __init__(self):
        """Build default object."""
        self.value = None
        self.lower_bound = None
        self.upper_bound = None
        self.variable_references = ListOfVariableReferences()
        self.parameter_references = ListOfParameterReferences()

    def is_empty(self):
        """Return whether all attributes are unspecified."""
        return (self.value is None and self.lower_bound is None and self.upper_bound is None) or self.variable_references.is_empty()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        if self.value is not None:
            result.set('value', str(self.value))
        if self.lower_bound is not None:
            result.set('lowerBound', str(self.lower_bound))
        if self.upper_bound is not None:
            result.set('upperBound', str(self.upper_bound))
        opt = [self.variable_references, self.parameter_references]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls()
        result._init_from_xml_node(node)
        n = get_unique_child(node, ListOfVariableReferences.tag, False)
        if n is not None:
            result.variable_references = ListOfVariableReferences.from_xml_node(n)
        n = get_unique_child(node, ListOfParameterReferences.tag, False)
        if n is not None:
            result.parameter_references = ListOfParameterReferences.from_xml_node(n)
        return result

    def _init_from_xml_node(self, node):
        """Match attributes with given node."""
        self.value = node.get('value')
        self.lower_bound = node.get('lowerBound')
        self.upper_bound = node.get('upperBound')


class VariableReference(object):
    """
    Reference to RBA-problem decision-variable in linear constraint

    Attributes
    ----------
    variable : str
        ID of decision variable in RBA-problem, in linear sum on LHS of constraint
    coefficient : float
        Linear coefficient, in linear sum on LHS of constraint
    """

    tag = 'variableReference'
    
    def __init__(self, variable, constant_coefficient, parameter_coefficient=None, multiplied_with_growth_rate=None):
        """
        Constructor.

        Parameters
        ----------
        variable : str
            ID of decision variable in RBA-problem, in linear sum on LHS of constraint
        coefficient : float
            Linear coefficient, in linear sum on LHS of constraint
        parameter_coefficient : str
            ID of parameter in RBA- parameters.xml file
        multiplied_with_growth_rate : bool
            True if linear coefficients should be proportional to (multiplied with) growth rate (False by default).
        """

        self.variable = variable
        self.constant_coefficient = constant_coefficient
        self.parameter_coefficient = parameter_coefficient
        self.multiplied_with_growth_rate = is_true(multiplied_with_growth_rate) if multiplied_with_growth_rate is not None else False

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('variable', self.variable)
        result.set('constant_coefficient', str(self.constant_coefficient))
        if self.parameter_coefficient is not None:
            result.set('parameter_coefficient', self.parameter_coefficient)
        result.set('multiplied_with_growth_rate', str(self.multiplied_with_growth_rate).lower())
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('variable'), float(node.get('constant_coefficient')),node.get('parameter_coefficient'),node.get('multiplied_with_growth_rate'))


class ParameterReference(object):
    """
    Reference to RBA model-parameter (µ-dependent, not decision variable dependent) in linear constraint

    Attributes
    ----------
    parameter : str
        ID of parameter in RBA- parameters.xml file
    coefficient : float
        Linear coefficient
    """

    tag = 'parameterReference'

    def __init__(self, parameter,constant_coefficient, multiplied_with_growth_rate=None):
        """
        Constructor.

        Parameters
        ----------
        species : str
            Identifier of species.
        stoichiometry : float
            Stoichiometry of species.
        multiplied_with_growth_rate : bool
            True if linear coefficients should be proportional to (multiplied with) growth rate (False by default).

        """
        self.parameter = parameter
        self.constant_coefficient = constant_coefficient
        self.multiplied_with_growth_rate = is_true(multiplied_with_growth_rate) if multiplied_with_growth_rate is not None else False

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('parameter', self.parameter)
        result.set('constant_coefficient', str(self.constant_coefficient))
        result.set('multiplied_with_growth_rate', str(self.multiplied_with_growth_rate).lower())
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('parameter'), float(node.get('constant_coefficient')), node.get('multiplied_with_growth_rate'))


class ConstraintReference(object):
    """
    Reference to RBA-problem constraint

    Attributes
    ----------
    constraint : str
        ID of constraint (row) in RBA-problem
    """

    tag = 'constraintReference'
    
    def __init__(self, constraint):
        """
        Constructor.

        Parameters
        ----------
        constraint : str
            ID of constraint (row) in RBA-problem
        """
        self.constraint = constraint

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('constraint', self.constraint)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('constraint'))


class ListOfConstraints(ListOf):
    """List of Constraint elements."""

    tag = 'listOfConstraints'
    list_element = Constraint


class ListOfVariableReferences(ListOf):
    """List of VariableReference elements."""

    tag = 'listOfVariableReferences'
    list_element = VariableReference


class ListOfConstraintReferences(ListOf):
    """List of ConstraintReference elements."""

    tag = 'listOfConstraintReferences'
    list_element = ConstraintReference


class ListOfParameterReferences(ListOf):
    """List of ParameterReference elements."""

    tag = 'listOfParameterReferences'
    list_element = ParameterReference

