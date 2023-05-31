import os
import argparse
import write_pdb


def main():
    """Check which transformations are needed, load, exec and save."""
    parser = argparse.ArgumentParser(
        description="Given a pdb file, applies requested transformations."
    )
    parser.add_argument("--file", "-f", help="load name")
    parser.add_argument("--save", "-s", help="save name")
    parser.add_argument(
        "--fix_atom_numbering",
        "-a",
        help="Fix the atom numbering in the pdb file.",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--fix_resiude_numbering",
        "-r",
        help="Fix the residue numbering in the pdb file.",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--section_into_chains",
        "-c",
        help="Section into chains. If integer > 0 given, sections into chains with that many res.",
        type=int,
        default=0,
    )
    args = parser.parse_args()

    if not os.path.isfile(args.file):
        raise FileNotFoundError("Please give .pdb input.")
    if not args.save:
        args.save = args.file
    with open(args.file, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    if args.fix_atom_numbering:
        lines = write_pdb.fix_atom_numbering(lines)
    if args.fix_resiude_numbering:
        lines = write_pdb.fix_atom_numbering(lines)
    if args.section_into_chains:
        lines = write_pdb.section_into_chains(lines, args.section_into_chains)
    with open(args.save, "w", encoding="utf-8") as fh:
        fh.writelines(lines)


if __name__ == "__main__":
    main()
