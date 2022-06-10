import sys
from write_pdb.pdbline import PDBLINE


def main():
    options = {

    }
    filename = sys.argv[1]

    lines = open(filename)
    lines = fix_atom_numbering(lines)
    savename = input(
        f'Loaded file {filename}, please enter save name'
        'or empty string to overwrite. Submit with ENTER.'
        )
    if savename is None:
        print('Overwriting.')
    else:
        print(f'Writing to {savename}.')
    with open(filename, 'w') as fh:
        fh.writelines(lines)


if __name__ == '__main__':
    main()