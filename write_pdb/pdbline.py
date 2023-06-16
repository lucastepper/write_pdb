import numpy as np


class PDBLINE(dict):
    """
    Create a single PDB line and write its fields.
    """
    def __init__(self, input_line, **kwargs):
        super().__init__()
        self.input_line = input_line
        # keep track if the line is an ATOM or HETATM
        self.is_atom = True
        # should actually make field_templates private Standart pdb
        # definition, changed to make it zero based and compatible
        # with python indexing by subtracting one from the left index
        self.field_templates = {
            'type':        ([0, 6],     '“ATOM”',		                     'character'),
            'atomid':      ([6, 11],	'Atom serial number	right',	         'integer'),
            'atomname':    ([12, 16],	'Atom name	left',	                 'character'),
            'altlocid':    ([16, 17],   'Alternate location indicator',	     'character'),
            'resname':     ([17, 20],	'Residue name	right',	             'character'),
            'chainid':     ([21, 22],   'Chain identifier',		             'character'),
            'resid':       ([22, 26],	'Residue sequence number	right',	 'integer'),
            'code':        ([26, 27],   'Code for insertions of residues',	 'character'),
            'posx':        ([30, 38],	'X orthogonal A coordinate	right',  'real (8.3)'),
            'posy':        ([38, 46],	'Y orthogonal A coordinate	right',  'real (8.3)'),
            'posz':        ([46, 54],	'Z orthogonal A coordinate	right',  'real (8.3)'),
            'occupancy':   ([54, 60],	'Occupancy	right',	                 'real (6.2)'),
            'tempfact':    ([60, 66],	'Temperature factor	right',	         'real (6.2)'),
            'segid':       ([72, 76],	'Segment identifier	left',	         'character'),
            'element':     ([76, 78],	'Element symbol	right',	             'character'),
        }
        for key, item in kwargs.items():
            self.__setitem__(key, item)

    @staticmethod
    def from_line(pdb_line):
        """ Construct a PDBLINE obj from a full pdb line given as a str. """
        # if not a proper atom line, mark and return
        if pdb_line[ : 4].strip() not in ['ATOM', 'HETATM']:
            # since I do not fix numbering on TER lines jet, just erase it
            if pdb_line.startswith('TER'):
                pdb_line = 'TER\n'
            to_construct = PDBLINE(pdb_line)
            to_construct.is_atom = False
            return to_construct
        # pad or cut line to len 79
        if len(pdb_line) < 79:
            pdb_line += ' ' * (79 - len(pdb_line))
        pdb_line = pdb_line[  : 79]
        to_construct = PDBLINE(pdb_line)
        for key, field_descriptor in to_construct.field_templates.items():
            idx_start, idx_end = field_descriptor[0]
            to_construct[key] = pdb_line[idx_start : idx_end]
        return to_construct

    def get_str_with_len(self, content, length):
        """ Take the content and pad it with spaces from the left until
        if has length. Asserts that str(content) <= length.
        """
        out_str = str(content)
        assert len(out_str) <= length, f'content too long: {content}, {length}'
        while len(out_str) < length:
            out_str = ' ' + out_str
        return out_str

    def __setitem__(self, key, item):
        idxs = self.field_templates[key][0]
        len_field = idxs[1] - idxs[0]
        assert len(str(item)) <= len_field, (
            'The given item to be written into a PDBLINE does not have the right '
            'length. When cast to a string, it should have less than or equal '
            f'than {len_field}, but it has: {len(str(item))}.'
        )
        super().__setitem__(key, self.get_str_with_len(item, len_field))

    def sub_line(self, idx1, idx2, word, line):
        """ Given a pdb line string, substitute a single fields with word. """
        assert idx2 - idx1 == len(word), f'{idx1 - idx2} {len(word)}'
        for i, idx in enumerate(range(idx1, idx2)):
            line[idx] = word[i]

    def get_line(self, pad_with=' \n'):
        """ Write all fields into a line, empty fields end up as spaces.
        If the line is not of type ATOM or HETATOM, return unmodified line. """
        if not self.is_atom:
            return self.input_line
        line = 79 * [' ']
        for key, (colums, __, __) in self.field_templates.items():
            if key in super().keys():
                word = self.get_str_with_len(super().__getitem__(key), colums[1] - colums[0])
                self.sub_line(colums[0], colums[1], word, line)
        return ''.join(line) + pad_with

    def set_positions(self, pos: np.ndarray):
        """ Set the position of the atom. """
        # I actually do not know what the pdb standart does with positions > 100
        self["posx"] = round(pos[0], 3)
        self["posy"] = round(pos[1], 3)
        self["posz"] = round(pos[2], 3)

    def __repr__(self):
        return self.get_line().strip('\n')
