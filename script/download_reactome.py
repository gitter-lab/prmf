#!/usr/bin/env python
import xml.etree.ElementTree as ET
import argparse
import itertools as it
import os, os.path
import requests

EXPRESSION_PATHWAYS_URL = "https://reactome.org/ContentService/data/pathway/R-HSA-74160/containedEvents"

# TODO can this be read from the xml?
BIOPAX_NS = '{http://www.biopax.org/release/biopax-level3.owl#}'
CONTROL = BIOPAX_NS + "Control"
CONTROL_TYPE = BIOPAX_NS + "controlType"
CONTROLLER = BIOPAX_NS + "controller"
CONTROLLED = BIOPAX_NS + "controlled"
COMPONENT = BIOPAX_NS + "component"

# entityReference and its subclasses
ENTITY_REFERENCE = BIOPAX_NS + "entityReference"
DNA_REFERENCE = BIOPAX_NS + "DnaReference"
PROTEIN_REFERENCE = BIOPAX_NS + "ProteinReference"
RNA_REFERENCE = BIOPAX_NS + "RnaReference"
SMALL_MOLECULE_REFERENCE = BIOPAX_NS + "SmallMoleculeReference"
ENTITY_REFERENCE_CLASSES = [ENTITY_REFERENCE, DNA_REFERENCE, PROTEIN_REFERENCE, RNA_REFERENCE, SMALL_MOLECULE_REFERENCE]

XREF = BIOPAX_NS + "xref"
PROTEIN = BIOPAX_NS + "Protein"
COMPLEX = BIOPAX_NS + "Complex"
BIOCHEMICAL_REACTION = BIOPAX_NS + "BiochemicalReaction"
LEFT = BIOPAX_NS + "left"
RIGHT = BIOPAX_NS + "right"
DB = BIOPAX_NS + "db"
BIOPAX_ID = BIOPAX_NS + "id"

RDF_NS = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
ID = RDF_NS + "ID"

class BioPaxDoc:
  def __init__(self, fp):
    self.fp = fp
    self.tree = ET.parse(fp)
    self.root = self.tree.getroot()

    # top-level elements in the BioPax namespace have an ID attribute that is frequently referenced
    self.id_to_elem = {}
    for elem in self.root.findall('.//*'):
      if BIOPAX_NS in elem.tag:
        id_v = None
        for k,v in elem.attrib.items():
          if k == ID:
            id_v = v
            break
        self.id_to_elem[id_v] = elem

  def parse_control_as_sif(self, control_elem):
    sif_strs = []
    edge_str = None
    control_type = self.get_control_type(control_elem)
    if control_type == "ACTIVATION":
      edge_str = "->"
    elif control_type == "INHIBITION":
      edge_str = "-|"

    if edge_str is not None:
      controller = self.get_controller(control_elem)
      # TODO what if controller is not a complex: use a generic higher level function to figure out what kind of element then dispatch to get_complex_protein_ids, get_biochem_rxn_protein_ids, etc.
      controller_proteins = self.get_complex_protein_ids(controller)
      controlled = self.get_controlled(control_elem)
      controlled_proteins = self.get_biochem_rxn_protein_ids(controlled)

      for controller_protein, controlled_protein in it.product(controller_proteins, controlled_proteins):
        sif_str = " ".join([controller_protein, edge_str, controlled_protein])
        sif_strs.append(sif_str)

    return sif_strs

  def get_control_type(self, control_elem):
    control_type = None
    control_type_elem = control_elem.find(CONTROL_TYPE)
    if control_type_elem is not None:
      control_type = control_type_elem.text
    return control_type

  # TODO biopax documentation says there can be 0 or more
  def get_control_member(self, control_elem, member_const):
    """
    Parameters
    ----------
    member_const : str
      one of CONTROLLER or CONTROLLED
    """
    member = None
    member_elem = control_elem.find(member_const)
    if member_elem is not None:
      member_id = None
      member = self.resolve_link_by_id(member_elem)
    return member

  def get_controller(self, control_elem):
    return self.get_control_member(control_elem, CONTROLLER)

  def get_controlled(self, control_elem):
    return self.get_control_member(control_elem, CONTROLLED)

  def get_complex_proteins(self, complex_elem):
    """
    Return all PROTEIN elements that are members of a complex
    """
    protein_elems = []
    complex_elems = []
    complex_elems.append(complex_elem)
    while len(complex_elems) != 0:
      complex_elem = complex_elems.pop()
      for component_link_elem in complex_elem.findall(COMPONENT):
        component_elem = self.resolve_link_by_id(component_link_elem)
        if component_elem.tag == COMPLEX:
          complex_elems.append(component_elem)
        elif component_elem.tag == PROTEIN:
          protein_elems.append(component_elem)
    return protein_elems

  def get_biochem_rxn_protein_ids(self, rxn_elem):
    """
    In BioPax, gene expression is represented as a BIOCHEMICAL_REACTION from a DNA element to a PROTEIN element
    """
    biochem_rxn_proteins = []
    # TODO check conversion direction?
    # require left is dna and right is protein
    left_elem_db_to_id = self.resolve_entity_reference(self.resolve_link_by_id(rxn_elem.find(LEFT)))
    right_elem_db_to_id = self.resolve_entity_reference(self.resolve_link_by_id(rxn_elem.find(RIGHT)))
    for db, ids in right_elem_db_to_id.items():
      for id_v in ids:
        biochem_rxn_proteins.append(id_v)
    return biochem_rxn_proteins

  def resolve_entity_reference(self, elem):
    """
    Return a mapping of database (namespace) to identifiers for an element with an ENTITY_REFERENCE or one of its subclasses
    """
    all_db_to_ids = {}
    for entity_reference in ENTITY_REFERENCE_CLASSES:
      for entity_ref_link_elem in elem.findall(entity_reference):
        entity_ref_elem = self.resolve_link_by_id(entity_ref_link_elem)
        db_to_id = self.resolve_xrefs(entity_ref_elem)
        for k,v in db_to_id.items():
          if k in all_db_to_ids:
            all_db_to_ids[k].append(v)
          else:
            all_db_to_ids[k] = [v]
    return all_db_to_ids

  def get_complex_protein_ids(self, complex_elem):
    """
    With each protein represented as a member in the returned list, return a mapping of database (namespace) to identifiers for all PROTEIN elements in a COMPLEX element
    """
    # TODO return protein_ids or protein_db_to_ids?
    all_protein_db_to_ids = []
    protein_ids = []
    protein_elems = self.get_complex_proteins(complex_elem)
    for protein_elem in protein_elems:
      protein_db_to_ids = self.resolve_entity_reference(protein_elem)
      for db, ids in protein_db_to_ids.items():
        for id_v in ids:
          protein_ids.append(id_v)
      all_protein_db_to_ids.append(protein_db_to_ids)
    return protein_ids

  def resolve_link_by_id(self, elem):
    linked_elem = None
    id_v = None
    # TODO what if multiple links
    for k,v in elem.attrib.items():
      if "#" in v:
        id_v = v.replace("#", "")
        break
    if id_v is not None:
      linked_elem = self.id_to_elem.get(id_v)
    return linked_elem

  def resolve_xrefs(self, elem):
    """
    Return a mapping of database to identifier value
    """
    db_to_id = {}
    for xref_link_elem in elem.findall(XREF):
      xref_elem = self.resolve_link_by_id(xref_link_elem)

      db = None
      db_elem = xref_elem.find(DB)
      if db_elem is not None:
        db = db_elem.text

      id_v = None
      id_elem = xref_elem.find(BIOPAX_ID)
      if id_elem is not None:
        id_v = id_elem.text

      db_to_id[db] = id_v

    return db_to_id

def get_expression_pathways(outdir):
  ofp = os.path.join(args.outdir, 'pathways.json')
  req_obj = requests.get(EXPRESSION_PATHWAYS_URL)

if __name__ == "__main__":
  # TODO update CLI
  parser = argparse.ArgumentParser()
  parser.add_argument("--infile", "-i", required=True)
  parser.add_argument("--outdir", "-o", required=True)
  args = parser.parse_args()

  bioPaxDoc = BioPaxDoc(args.infile)
  for control_elem in bioPaxDoc.root.findall('.//' + CONTROL):
    for sif_str in bioPaxDoc.parse_control_as_sif(control_elem):
      print(sif_str)

