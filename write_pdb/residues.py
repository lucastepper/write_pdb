# pylint: disable=line-too-long
from . import Residue


# to do: add the rest of the canonical aas
ACE = Residue('ACE', 'ATOM',
    ['1HH3', 'CH3', '2HH3', '3HH3', 'C', 'O'],
    ['H', 'C', 'H', 'H', 'C', 'O'])
NME = Residue('NME', 'ATOM',
    ['N', 'H', 'CH3', '1HH3', '2HH3', '3HH3'],
    ['N', 'H', 'C', 'H', 'H', 'H'])
ALA = Residue('ALA','ATOM',
    ['N', 'H', 'CA', 'HA', 'CB', '1HB', '2HB', '3HB', 'C', 'O'],
    ['N', 'H', 'C', 'H', 'C', 'H', 'H', 'H', 'C', 'O'])
ASN = Residue('ASN', 'ATOM',
    ['N', 'H1', 'H2', 'H3', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O'],
    ['N', 'H', 'H', 'H', 'C', 'H', 'C', 'H', 'H', 'C', 'O', 'N', 'H', 'H', 'C', 'O'])
PHE = Residue('PHE', 'ATOM',
    ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O'],
    ['N','H','C','H','C','H','H','C','C','H','C','H','C','H','C','H','C','H','C','O'])
GLY = Residue('GLY', 'ATOM',
    ['N','H','CA','HA1','HA2','C','O'],
    ['N', 'H', 'C', 'H', 'H', 'C', 'O'])
ILE = Residue('ILE', 'ATOM',
    ['N', 'H', 'CA', 'HA', 'CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG11', 'HG12', 'CD', 'HD1', 'HD2', 'HD3', 'C', 'O'],
    ['N', 'N', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'H', 'C', 'O'])
SER = Residue('SER', 'ATOM',
    ['N', 'H', 'C', 'H', 'C', 'H', 'H', 'O', 'H', 'C', 'O'],
    ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'OG', 'HG', 'C', 'O'])
LEU = Residue('LEU', 'ATOM',
    ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O'],
    ['N' , 'H' , 'C' , 'H' , 'C' , 'H' , 'H' , 'C' , 'H' , 'C' , 'H' , 'H' , 'H' , 'C' , 'H' , 'H' , 'H' , 'C' , 'O'])
