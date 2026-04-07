"""Module defining density-specific classes for RBA XML structures."""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml._rbaml_version import __rbaml_version__ as rbaml_version
from rba.xml.common import get_unique_child, ListOf, is_true, xml_input_tag_error

__all__ = ['RbaCompartments', 'Compartment','CompartmentComposition','ConstituentReference', 'ListOfCompartments','ListOfConstituentReferences']

class RbaCompartments(object):
    """
    Compartment-related structures.

    Attributes
    ----------
    compartments : ListOfCompartments
        List of cell compartments.
    constraints : ListOfLinearConstraints
        List of user defined custom constraints to add.
    """
    tag='RBACompartments'
    def __init__(self):
        """Constructor."""
        self.compartments = ListOfCompartments()

    def write(self, output_stream,meta_data={}):
        """
        Write information as an XML structure.

        Parameters
        ----------
        output_stream : file or buffer
            Location where XML structure should be written.
        """
        root = etree.Element(self.tag)
        root.extend([self.compartments.to_xml_node()])
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
        # Verify that root tag in file is the one specified by class RBACompartments
        if node.tag != cls.tag:
            raise xml_input_tag_error(file=input_stream.name,file_tag=node.tag,requested_tag=cls.tag)

        result = cls()
        n = get_unique_child(node, ListOfCompartments.tag)
        result.compartments = ListOfCompartments.from_xml_node(n)
        return result


class Compartment(object):
    """
    Compartment information.

    Attributes
    ----------
    id : str
        identifier of compartment.
    """

    tag = 'compartment'

    def __init__(self, id_,lb_=None,ub_=None,is_external_=False):
        """
        Constructor.

        Parameters
        ----------
        id : str
            identifier of compartment.
        """
        self.id = id_
        self.lower_bound = lb_
        self.upper_bound = ub_
        self.is_external = is_external_ if is_external_ is not None else False
        self.composition = CompartmentComposition()

    def to_xml_node(self):
        """
        Convert to xml node
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        if self.lower_bound is not None:
            result.set('lowerBound', str(self.lower_bound))
        if self.upper_bound is not None:
            result.set('upperBound', str(self.upper_bound))
        result.set('isExternal', str(self.is_external).lower())
        if not self.composition.is_empty():
            result.extend([self.composition.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result =  cls(node.get('id'),node.get('lowerBound'),node.get('upperBound'),is_true(node.get('isExternal')))
        n = get_unique_child(node, CompartmentComposition.tag, False)
        if n is not None:
            result.composition = CompartmentComposition.from_xml_node(n)
        return result


class CompartmentComposition(object):
    """
    compartment composition information.

    Attributes
    ----------
    """

    tag = 'compartmentComposition'

    def __init__(self):
        """
        Constructor.

        Parameters
        ----------
        id : str
            identifier of compartment.
        """
        self.constituents = ListOfConstituentReferences()

    def is_empty(self):
        """Return whether composition is fully empty."""
        return self.constituents.is_empty()

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        # optional subelements
        opt = [self.constituents]
        result.extend([e.to_xml_node() for e in opt if not e.is_empty()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        result = cls()
        n = get_unique_child(node, ListOfConstituentReferences.tag, False)
        if n is not None:
            result.constituents = ListOfConstituentReferences.from_xml_node(n)
        return result


class ConstituentReference(object):
    """
    Reference compartment constituent / housekeeping species (macromolecule)

    Attributes
    ----------
    species : str
        ID of macromolecule in RBA-problem
    fraction : str
        occupation/weight fraction of macromolecule in compartment
    """

    tag = 'constituentReference'
    
    def __init__(self, species,fraction=None):
        """
        Constructor.

        Parameters
        ----------
        constraint : str
            ID of constraint (row) in RBA-problem
        """
        self.species = species
        self.fraction = fraction

    def to_xml_node(self):
        """Convert to xml node."""
        result = etree.Element(self.tag)
        result.set('species', self.species)
        if self.fraction is not None:
            result.set('fraction', self.fraction)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """Build object from xml node."""
        return cls(node.get('species'),node.get('fraction'))


class ListOfCompartments(ListOf):
    """
    List of Compartment elements
    """

    tag = 'listOfCompartments'
    list_element = Compartment


class ListOfConstituentReferences(ListOf):
    """
    List of ConstituentReference elements
    """

    tag = 'listOfConstituentReferences'
    list_element = ConstituentReference
