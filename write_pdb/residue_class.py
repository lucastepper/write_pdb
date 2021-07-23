from write_pdb.pdbline import PDBLINE


class Residue:
    """
    A class representing a aminoacid residue in a polypeptide.
    """
    def __init__(self, resname, atom_categroy, atomnames, elements):
        """ Constructor for a specific residue, ie alanine.
        Arguments:
            resname(str): Name for the residue.
            atom_categroy (str): Category for the atoms in the residue, ie ATOM or HETATOM.
            atom_names (iterable): Names for the atoms in the residue.
            elements (iterable): Element names for the atoms in the residue.
            atom_counter (int): Counter for the atom number in the entire chain, which
                the residue belongs to.
        """
        self.atom_categroy = atom_categroy
        self.atomnames = atomnames
        self.resname = resname
        self.elements = elements
        temp = PDBLINE()
        self.fields = temp.field_templates
    
    def sub_line(self, idx1, idx2, word, line):
        assert idx2 - idx1 + 1 == len(word), f'{idx1 - idx2 + 1} {len(word)}'
        for i, idx in enumerate(range(idx1, idx2 + 1)):
            line[idx - 1] = word[i]
    
    def get_str_with_len(self, content, length):
        out_str = str(content)
        assert len(out_str) <= length, f'content too long: {content}, {length}'
        while len(out_str) < length:
            out_str = ' ' + out_str
        return out_str
    
    def get_line(self, to_add):
        line = 79 * [' ']
        for key, (colums, __, __) in self.fields.items():
            if key in to_add:
                word = self.get_str_with_len(to_add[key], colums[1] - colums[0] + 1)
                self.sub_line(colums[0], colums[1], word, line)
        return ''.join(line) + '\n'
    
    def get_lines(self, resid, atomid, positions, chainid='A', first_res=False, last_res=False):
        lines = []
        counter = 0
        for i in range(len(self.atom_categroy)):
            try:
                lines.append(self.get_line({
                    'type': self.atom_categroy,
                    'atomid': counter + atomid,
                    'atomname': self.atomnames[i],
                    'resid': resid,
                    'resname': self.resname,
                    'chainid': chainid,
                    'posx': f'{positions[counter, 0]:.3f}',
                    'posy': f'{positions[counter, 1]:.3f}',
                    'posz': f'{positions[counter, 2]:.3f}',
                    'element': self.elements[i],
                }))
            except IndexError:
                print(positions, counter, positions.shape)
            counter += 1
            # add charge to Nitrogen, if first residue
            if first_res and i == 2:
                lines.append(self.get_line({
                    'type': self.atom_categroy,
                    'atomid': counter + atomid,
                    'atomname': 'H3',
                    'resid': resid,
                    'resname': self.resname,
                    'chainid': chainid,
                    'posx': f'{positions[counter, 0]:.3f}',
                    'posy': f'{positions[counter, 1]:.3f}',
                    'posz': f'{positions[counter, 2]:.3f}',
                    'element': 'H',
                }))
                counter += 1
        # add charge to Oxigen, if last residue
        if last_res:
            counter -= 1
            lines.pop(-1)
            for j in range(2):
                lines.append(self.get_line({
                    'type': self.atom_categroy,
                    'atomid': counter + atomid,
                    'atomname': f'OC{j}',
                    'resid': resid,
                    'resname': self.resname,
                    'chainid': chainid,
                    'posx': f'{positions[counter, 0]:.3f}',
                    'posy': f'{positions[counter, 1]:.3f}',
                    'posz': f'{positions[counter, 2]:.3f}',
                    'element': 'O',
                }))
                counter += 1
        return lines
