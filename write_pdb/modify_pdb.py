import string
from write_pdb.pdbline import PDBLINE


class IsNewChain():
    """ Small class to check if given PDBLine object is the
    same or different chain to the last one passed, based on chainid
    If not PDBLine.is_atom, returns False
    If the first PDBLine that is_atom is passed, returns return_val_first_call. """

    def __init__(self, return_val_first_call=True):
        self.chainid = None
        self.return_val_first_call = return_val_first_call

    def clear(self):
        self.chainid = None

    def __call__(self, pdb_line):
        if not pdb_line.is_atom:
            return False
        if self.chainid != pdb_line['chainid']:
            self.chainid = pdb_line['chainid']
            return self.return_val_first_call
        return False


def fix_atom_numbering(pdb_lines, restart_atomid_per_chain=True):
    """ Check that the atom numbering is increasing by one each line,
    given a list of pdb lines. Returns a bopy of pdb_lines.
    Arguments:
        pdblines (list of str): List containing the pdb lines to be changed.
        restart_atomid_per_chain (bool): Restart the atom numbering at 1 when a
            new chain begins or number starting at the beginning of the file to the end
            default: True
    """
    output_lines = []
    atom_counter = 1
    is_new_chain = IsNewChain(return_val_first_call=False)
    for line in pdb_lines:
        line_obj = PDBLINE.from_line(line)
        if restart_atomid_per_chain and is_new_chain(line_obj):
            atom_counter = 1
        if line_obj.is_atom:
            line_obj['atomid'] = atom_counter
            atom_counter += 1
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
    # store the name of resiude we are iterating over right now
    # to detect when it changes
    resname_cache = None
    for line in pdb_lines:
        line_obj = PDBLINE.from_line(line)
        # ignore REMARK TER and other such lines
        if not line_obj.is_atom:
            continue
        # set cache the first time a line belonging to residue is encountered
        if resname_cache is None:
            resname_cache = PDBLINE.from_line(line)['resname']
            continue
        # check if a new residue has started, count the number of residues
        if line_obj['resname'] != resname_cache:
            res_counter += 1
            resname_cache = line_obj['resname']
        line_obj['resid'] = res_counter
        output_lines.append(line_obj.get_line())
    return output_lines


def write_positions(pdb_lines, positions, idx_start=0, idx_end='inf'):
    """ Write new positions into a pdb file.
    Arguments
        pdblines (list of str): List containing the pdb lines to be changed.
        positions (np.ndarray (-1, 3)): Positions to be written into the pdb.
        idx_start (int): Index where to start changing pos in lines, 0-based
            default: 0; incluse
        idx_end (int): Index where to end changing pos in lines, 0-based
            default: int(1e99); exclusive
    """
    if idx_end == 'inf':
        idx_end = len(lines)
    output_lines = pdb_lines[ : idx_start]
    assert (idx_end - idx_start, 3) == positions.shape
    for line, pos in zip(pdb_lines[idx_start : idx_end], positions):
        line_obj = PDBLINE.from_line(line)
        # I actually do not know what the pdb standart does with positions > 100
        line_obj['posx'] = round(pos[0], 3)
        line_obj['posy'] = round(pos[1], 3)
        line_obj['posz'] = round(pos[2], 3)
        output_lines.append(line_obj.get_line())
    output_lines.extend(pdb_lines[idx_end : ])
    return output_lines


class IsNewResidue():
    """ Small class to check if given PDBLine object is the
    same or different residue to the last one passed.
    Res is different when either resid or resname changes
    If not PDBLine.is_atom, returns False
    If the first PDBLine that is_atom is passed, returns True. """

    def __init__(self):
        """ Init in empty state, set the first pdb_line input as ref. """
        self.resname = None
        self.resid = None

    def clear(self):
        """ Clear the ref, which we do ie. to start a new residue. """
        self.resname = None
        self.resid = None

    def __call__(self, pdb_line):
        if not pdb_line.is_atom:
            return False
        if self.resid != pdb_line['resid'] or self.resid != pdb_line['resid']:
            self.resname = pdb_line['resname']
            self.resid = pdb_line['resid']
            return True
        return False

def section_into_chains(pdb_lines, residues_per_chain, chain_names=string.ascii_uppercase):
    """ Define a given number of residues into a chain. Chains are identified
    by upper case letters in alphabetical order. Does not change residue numbering.
    Checks that a 'TER' line is inserted after a chain finishes.
    Arguments
        residues_per_chain (list of ints or int): Number of residues per chain.
            If int, all chains get the same number of residues
            If list, runs through the list and assigns a chain for each
            integer in the list with that number of residues.
        chain_names (list of str): List of chain names to be used. Must be at least
            as long as the number of chains that are defined by n_lines / residues_per_chain.
            default: string.ascii_uppercase
    """
    # we count the number of residues by iterating through all lines
    # where a new residue starts when either the resid or the resname changes
    is_new_residue = IsNewResidue()
    n_residues = int(sum(is_new_residue(PDBLINE.from_line(line)) for line in pdb_lines))
    print(f"Input {residues_per_chain=}, computed {n_residues=}")
    if n_residues > len(chain_names):
        raise ValueError("We have more residues than chain names, please provide more chain names")
    is_new_residue.clear()
    # set up residues per chain
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
    # now do the actual sectioning into chains
    output_lines = []
    chain_counter = 0
    res_counter = 0
    for i, line in enumerate(pdb_lines):
        line_obj = PDBLINE.from_line(line)
        # erase counter in TER lines
        if line.startswith('TER'):
            output_lines.append('TER \n')
        # leave other non ATOM HEATATM lines unchanged
        elif not line_obj.is_atom:
            output_lines.append(line)
        else:
            res_counter += is_new_residue(line_obj)
            # check if a new chain should start
            if int(res_counter) > sum(residues_per_chain[ : chain_counter + 1]):
                chain_counter += 1
                # check if the next line is a "TER" line, if not, add one
                if not "TER" in pdb_lines[i + 1]:
                    output_lines.append("TER \n")
            line_obj['chainid'] = chain_names[chain_counter]
            output_lines.append(line_obj.get_line())
    return output_lines
