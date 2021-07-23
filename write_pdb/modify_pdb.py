import string
from write_pdb.pdbline import PDBLINE


def fix_atom_numbering(pdb_lines):
    """ Check that the atom numbering is increasing by one each line,
    given a list of pdb lines. Returns a bopy of pdb_lines.
    Arguments:
        pdblines (list of str): List containing the pdb lines to be changed.
    """
    output_lines = []
    for i, line in enumerate(pdb_lines):
        line_obj = PDBLINE.from_line(line)
        line_obj['atomid'] = i
        output_lines.append(line_obj.get_line())
    return output_lines


def fix_residue_numbering(pdb_lines):
    """ Check that the residue numbering is increasing by one when the resname
    field chainges, given a list of pdb lines. Returns a bopy of pdb_lines.
    Arguments:
        pdblines (list of str): List containing the pdb lines to be changed.
    """
    output_lines = []
    res_counter = 1
    resname_cache = PDBLINE.from_line(pdb_lines[0])['resname']
    for line in pdb_lines:
        line_obj = PDBLINE.from_line(line)
        # check if a new residue has started, count the number of residues
        if line_obj['resname'] != resname_cache:
            res_counter += 1
            resname_cache = line_obj['resname']
        line_obj['resid'] = res_counter
        output_lines.append(line_obj.get_line())
    return output_lines


def change_positions(pdb_lines, positions, idx_start=0, idx_end=int(1e99)):
    """ Write new positions into a pdb file.
    Arguments
        pdblines (list of str): List containing the pdb lines to be changed.
        positions (np.ndarray (-1, 3)): Positions to be written into the pdb.
        idx_start (int): Index where to start changing pos in lines, 0-based
            default: 0; incluse
        idx_end (int): Index where to end changing pos in lines, 0-based
            default: int(1e99); exclusive
    """
    output_lines = pdb_lines[ : idx_start]
    assert (idx_end - idx_start - 1, 3) == positions.shape
    for line, pos in zip(pdb_lines[idx_start : idx_end], positions):
        line_obj = PDBLINE.from_line(line)
        # I actually do not know what the pdb standart does with positions > 100
        line_obj['posx'] = round(pos[0], 3)
        line_obj['posy'] = round(pos[1], 3)
        line_obj['posz'] = round(pos[2], 3)
        output_lines.append(line_obj.get_line())
    output_lines.extend(pdb_lines[idx_end : ])
    return output_lines


def section_into_chains(pdb_lines, residues_per_chain):
    """ Define a given number of residues into a chain. Chains are identified
    by upper case letters in alphabetical order. fixes the residue numbering
    Arguments
        residues_per_chain (list of ints or int): Number of residues per chain.
        If int, all chains get the same number of residues
        If list, runs through the list and assigns a chain for each integer in
        the list with that number of residues.
    """
    temp_lines = fix_residue_numbering(pdb_lines)
    n_residues = int(PDBLINE.from_line(temp_lines[-1])['resid'])
    if isinstance(residues_per_chain, int):
        assert n_residues % residues_per_chain == 0, (
            'The number of residues per chain is an integer, so every chain should '
            'have the same number of residues. residues_per_chain was given as '
            f'{residues_per_chain}, which is not a divisor of the number of residues'
            f'in pdb lines, which is {n_residues}. '
            f'{n_residues % residues_per_chain} are left hanging. '
        )
        residues_per_chain = n_residues // residues_per_chain * [residues_per_chain]
    elif isinstance(residues_per_chain, list):
        if sum([isinstance(x, int) for x in residues_per_chain]) != 0:
            raise TypeError('Type of residues_per_chain needs to be int or list of ints')
        n_residues_chain = sum(residues_per_chain)
        assert n_residues % pdb_lines == 0, (
            'The number of residues per chain is a list, so every chain should '
            'have the number of residues given in the list. The sum of integers '
            f'is f{n_residues_chain}, which is not the same as the total number of '
            f'in pdb_lines, which is {n_residues}'
        )
    else:
        raise TypeError('Type of residues_per_chain needs to be int or list of ints')

    output_lines = []
    chain_counter = 0
    for line in temp_lines:
        line_obj = PDBLINE.from_line(line)
        if int(line_obj['resid']) > sum(residues_per_chain[ : chain_counter + 1]):
            chain_counter += 1
        line_obj['chainid'] = string.ascii_uppercase[chain_counter]
        output_lines.append(line_obj.get_line())
    return output_lines
