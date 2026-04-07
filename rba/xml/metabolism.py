"""
Metabolism-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml._rbaml_version import __rbaml_version__ as rbaml_version
from rba.xml.common import (is_true, get_unique_child,
                            ListOf, ListOfProducts, 
                            ListOfReactants, xml_input_tag_error)

__all__ = ['RbaMetabolism', 
           'CompartmentMetabolism', 
           'ListOfCompartmentsMetabolism',
           'Species', 'ListOfSpecies', 'Reaction', 'ListOfReactions']


class RbaMetabolism(object):
    """
    Metabolism-related structures.

    Attributes
    ----------
    species : ListOfSpecies
        List of metabolites.
    reactions : ListOfReactions
        List of reactions.
    """
    tag='RBAMetabolism'
    def __init__(self):
        """
        Default constructor.
        """
        self.species = ListOfSpecies()
        self.reactions = ListOfReactions()

    def write(self, output_stream,meta_data={}):
        """
        Write information as an XML structure.

        Parameters
        ----------
        output_stream : file or buffer
            Location where XML structure should be written.
        """
        root = etree.Element(self.tag)
        root.extend([self.species.to_xml_node(), 
                     self.reactions.to_xml_node()])
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
        # Verify that root tag in file is the one specified by class RBAMetabolism
        if node.tag != cls.tag:
            raise xml_input_tag_error(file=input_stream.name,file_tag=node.tag,requested_tag=cls.tag)
        result = cls()
        try:
            #OBSOLETE! But still required to read and update models in old format (rbaml v1)
            n = get_unique_child(node, ListOfCompartmentsMetabolism.tag)
            result.compartments = ListOfCompartmentsMetabolism.from_xml_node(n)
        except: 
            pass
        n = get_unique_child(node, ListOfSpecies.tag)
        result.species = ListOfSpecies.from_xml_node(n)
        n = get_unique_child(node, ListOfReactions.tag)
        result.reactions = ListOfReactions.from_xml_node(n)
        return result
#

class CompartmentMetabolism(object):
    """
    Compartment information originally located in RBAmetabolism..
    OBSOLETE! But still required to read and update models in old format (rbaml v1)

    Attributes
    ----------
    id : str
        identifier of compartment.
    """

    tag = 'compartment'

    def __init__(self, id_):
        """
        Constructor.

        Parameters
        ----------
        id : str
            identifier of compartment.
        """
        self.id = id_

    def to_xml_node(self):
        """
        Convert to xml node
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('id'))


class ListOfCompartmentsMetabolism(ListOf):
    """
    List of Compartment elements originally located in RBAmetabolism.
    OBSOLETE! But still required to read and update models in old format (rbaml v1)
    """

    tag = 'listOfCompartments'
    list_element = CompartmentMetabolism


class Species(object):
    """
    Chemical species.

    Attributes
    ----------
    id : str
        Identifier of species.
    boundary_condition : bool
        Whether the species belongs to the boundary of the system.
    """

    tag = 'species'

    def __init__(self, id_, boundary_condition):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier of species.
        boundary_condition: bool
            Whether the species belongs to the boundary of the system.
        """
        self.id = id_
        self.boundary_condition = boundary_condition

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('boundaryCondition', str(self.boundary_condition).lower())
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('id'), is_true(node.get('boundaryCondition')))


class ListOfSpecies(ListOf):
    """
    List of Species elements.
    """

    tag = 'listOfSpecies'
    list_element = Species


class Reaction(object):
    """
    Reaction.

    Attributes
    ----------
    id : str
        Identifier.
    reversible : bool
        True if reaction is reversible.
    reactants : ListOfReactants
        List of chemicals consumed by reaction.
    products : ListOfProducts
        List of chemicals produced by reaction.
    """

    tag = 'reaction'

    def __init__(self, id_, reversible):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        reversible : bool
            True if reaction is reversible.
        """
        self.id = id_
        self.reversible = reversible
        self.reactants = ListOfReactants()
        self.products = ListOfProducts()

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('reversible', str(self.reversible).lower())
        result.extend([self.reactants.to_xml_node(),
                       self.products.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(node.get('id'), is_true(node.get('reversible')))
        n = get_unique_child(node, ListOfReactants.tag, False)
        if n is not None:
            result.reactants = ListOfReactants.from_xml_node(n)
        n = get_unique_child(node, ListOfProducts.tag, False)
        if n is not None:
            result.products = ListOfProducts.from_xml_node(n)
        return result


class ListOfReactions(ListOf):
    """
    List of Reaction elements.
    """

    tag = 'listOfReactions'
    list_element = Reaction
