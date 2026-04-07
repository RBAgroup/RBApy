"""
Module defining macromolecule-specific classes used for RBA XML structures.
"""

# python 2/3 compatiblity
from __future__ import division, print_function, absolute_import

# global imports
from lxml import etree

# local imports
from rba.xml._rbaml_version import __rbaml_version__ as rbaml_version
from rba.xml.common import get_unique_child, ListOf, xml_input_tag_error

__all__ = ['RbaProteins','RbaRNAs','RbaDNA','RbaMacromolecules', 'Component', 'ListOfComponents',
           'Macromolecule', 'ListOfMacromolecules', 'ComponentReference',
           'Composition']


class RbaMacromolecules(object):
    """
    Macromolecule-related structures.

    Attributes
    ----------
    components : ListOfComponents
        List of components of macromolecule (e.g. amino acids for proteins).
    macromolecules : ListOfMacromolecules
        List of macromolecules.
    """
    tag='RBAMacromolecules'
    def __init__(self):
        """
        Default constructor.
        """
        self.components = ListOfComponents()
        self.macromolecules = ListOfMacromolecules()

    def write(self, output_stream,meta_data={}):
        """
        Write information as an XML structure.

        Parameters
        ----------
        output_stream : file or buffer
            Location where XML structure should be written.
        """
        root = etree.Element(self.tag)
        root.extend([self.components.to_xml_node(),
                     self.macromolecules.to_xml_node()])
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
        # Verify that root tag in file is the one specified by class RBAMacromolecules,RBAProteins,RBARnas,RBADna
        if node.tag != cls.tag:
            raise xml_input_tag_error(file=input_stream.name,file_tag=node.tag,requested_tag=cls.tag)
        result = cls()
        n = get_unique_child(node, ListOfComponents.tag)
        result.components = ListOfComponents.from_xml_node(n)
        n = get_unique_child(node, ListOfMacromolecules.tag)
        result.macromolecules = ListOfMacromolecules.from_xml_node(n)
        return result


class Component(object):
    """
    Component of a macromolecule.

    Attributes
    ----------
    id : str
        Identifier.
    name : str
        Usual name.
    type : str
        Type of molecule (amino acid, vitamin, etc.)
    weight : float
        Weight.
    """

    tag = 'component'

    def __init__(self, id_, name, type_, weight):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        name : str
            Usual name.
        type_ : str
            Type of molecule (amino acid, vitamin, etc.)
        weight : float
            Weight.
        """
        self.id = id_
        self.name = name
        self.type = type_
        self.weight = weight

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('name', self.name)
        result.set('type', self.type)
        result.set('weight', str(self.weight))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('id'), node.get('name'),
                   node.get('type'), node.get('weight'))


class ListOfComponents(ListOf):
    """
    List of Component elements.
    """

    tag = 'listOfComponents'
    list_element = Component


class Macromolecule(object):
    """
    Macromolecule.

    Attributes
    ----------
    id : str
        Identifier.
    compartment : str
        Identifier of compartment where molecule lives.
    composition : Composition
        Composition of macromolecule in terms of components.
    half_life : str or None
        Reference to parameter, representing the half_life time 
        of respective macromolecule. 

    """

    tag = 'macromolecule'

    def __init__(self, id_, compartment, composition=None, half_life=None):
        """
        Constructor.

        Parameters
        ----------
        id_ : str
            Identifier.
        compartment : str
            Identifier of compartment where molecule lives.
        composition : dict, optional
            Dictionary where keys are ids of components and values are their
            stoichiometry within the molecule.
        half_life : str , optional
            Peference to parameter (parameter-ID as str), representing the half_life time 
            of respective macromolecule. Defeault: None (no decay assumed)

        """
        self.id = id_
        self.compartment = compartment
        self.half_life = half_life
        self.composition = Composition()
        if composition:
            for comp, sto in composition.items():
                self.composition.append(ComponentReference(comp, sto))

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('id', self.id)
        result.set('compartment', self.compartment)
        if self.half_life is not None:
            result.set('half_life', self.half_life)
        if not self.composition.is_empty():
            result.extend([self.composition.to_xml_node()])
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        result = cls(id_=node.get('id'), compartment=node.get('compartment'), half_life=node.get('half_life'))
        comp_node = get_unique_child(node, Composition.tag)
        result.composition = Composition.from_xml_node(comp_node)
        return result


class ListOfMacromolecules(ListOf):
    """
    List of Macromolecule elements.
    """

    tag = 'listOfMacromolecules'
    list_element = Macromolecule


class ComponentReference(object):
    """
    Reference to a component, including stoichiometry.

    Attributes
    ----------
    component : str
        Identifier of component.
    stoichiometry : stoichiometry
        Stoichiometry of component.
    """

    tag = 'componentReference'

    def __init__(self, component, stoichiometry):
        """
        Constructor.

        Parameters
        ----------
        component : str
            Identifier of component.
        stoichiometry : stoichiometry
            Stoichiometry of component.
        """
        self.component = component
        self.stoichiometry = stoichiometry

    def to_xml_node(self):
        """
        Convert to xml node.
        """
        result = etree.Element(self.tag)
        result.set('component', self.component)
        result.set('stoichiometry', str(self.stoichiometry))
        return result

    @classmethod
    def from_xml_node(cls, node):
        """
        Constructor from xml node.
        """
        return cls(node.get('component'), float(node.get('stoichiometry')))


class Composition(ListOf):
    """
    List of ComponentReference elements.
    """

    tag = 'composition'
    list_element = ComponentReference


class RbaProteins(RbaMacromolecules):
    tag='RBAProteins'

class RbaDNA(RbaMacromolecules):
    tag='RBADna'

class RbaRNAs(RbaMacromolecules):
    tag='RBARnas'
