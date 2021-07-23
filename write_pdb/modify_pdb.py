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
