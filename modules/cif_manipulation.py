from __future__ import annotations

import ccdc
from ccdc import io
from ccdc.crystal import Crystal

from multiprocessing import current_process as cpr

# from .log_utils import get_logger


def readentry(input_cif: str) -> Crystal:
    """
    Reads a CIF file containing molecular structure data and converts it to a
    standard atom labeling convention using the ccdc.crystal module.

    Parameters:
        input_cif (str): filename (.CIF) containing crystal structure data.

    Returns:
        newcif (ccdc.crystal.Crystal): Crystal object containing structural data
                                       in the standard atom labeling convention.
    """
    # read in the cif to a crystal object
    with io.CrystalReader(input_cif, format="cif") as readcif:
        cif = readcif[0]
        cif.assign_bonds()
    readcif.close()

    # to remove duplicate atoms, need the empirical formula
    formula = cif.formula
    elamnt = formula.split(" ")

    # now convert to standard labelling convention and identify
    # duplicate atoms to be removed
    with open(input_cif, "r") as file:
        file.seek(0)
        newstring = str()
        lines = file.readlines()
        loop_pos = 0
        start = 0
        end = 0
        columncount = 0

        for i, line in enumerate(lines):
            # locate atom type and site label columns
            if "loop_" in line:
                loop_pos = i
            if ("_atom" in line) and ("_geom" not in line) and ("_aniso" not in line):
                start = loop_pos + 1
                end = i + 1
        for i in range(start, end):
            if "atom_site_type_symbol" in lines[i]:
                type_pos = columncount
            if "atom_site_label" in lines[i]:
                label_pos = columncount
            columncount += 1
        # keep track of relative position of type and label
        counting = {}
        cutoff = {}
        # to_remove = []
        for i in range(end, len(lines)):
            if "loop_" in lines[i]:
                break
            # lines with atom information will contain a ., so only look at these
            if "." in lines[i]:
                # split lines by whitespace
                col = lines[i].split()
                # keep count of how many of each element type
                if col[type_pos] not in counting:
                    counting[col[type_pos]] = 1
                elif col[type_pos] in counting:
                    counting[col[type_pos]] += 1
                # new atom labels
                newlabel = f"{col[type_pos]}{counting[col[type_pos]]}"
                lines[i] = lines[i].replace(col[label_pos], newlabel)
                # cutoff repeated atoms
                if newlabel in elamnt:
                    cutoff[col[type_pos]] = counting[col[type_pos]]

        # combine to new string
        for i in lines:
            newstring += i
        # read into new crystal object and assign bonds
        newcif = Crystal.from_string(newstring, format="cif")
        newcif.assign_bonds()
        file.close()
    return newcif


def get_coordinates(
    molecule: ccdc.molecule.Molecule, solvents: list[str]
) -> list[list[float]]:
    """
    Gets the [x,y,z] coordinate lists of the individual atoms due to be removed.

    Parameters:
        molecule (ccdc.molecule.Molecule): Molecule object representing full structure.
        solvents (list[str]): list of atom labels of solvents to be removed.

    Returns:
        atoms_coordinates (list[list[float]]): list of lists containing the 3 fractional
                                               coordinates [x, y, z] of each atom due for
                                               removal.

    """
    atoms_coordinates = []
    for atom in molecule.atoms:
        if atom.label in solvents:
            coord_unit = []
            coord_unit.append(atom.fractional_coordinates.x)
            coord_unit.append(atom.fractional_coordinates.y)
            coord_unit.append(atom.fractional_coordinates.z)
            # firmatting the coordinates to match the cif format
            for i in range(len(coord_unit)):
                if 0 < coord_unit[i] < 1:
                    coord_unit[i] = format(coord_unit[i], ".5f")
                elif coord_unit[i] == -0.0:
                    coord_unit[i] = "0.00000"
                elif coord_unit[i] < 0:
                    if coord_unit[i] >= -1:
                        coord_unit[i] += 1
                    elif coord_unit[i] < -1:
                        coord_unit[i] += 2
                    coord_unit[i] = format(coord_unit[i], ".5f")
                elif coord_unit[i] >= 1:
                    if coord_unit[i] >= 2:
                        coord_unit[i] -= 2
                    else:
                        coord_unit[i] -= 1
                    coord_unit[i] = format(coord_unit[i], ".5f")
            # writing the individual coordinate list to the total list of coordinates
            atoms_coordinates.append(coord_unit)
    return atoms_coordinates


def remove_solvents_from_file(
    lines: list[str], coordinates: list[list[float]]
) -> tuple[list[str], int]:
    """
    Parses CIF file and removes lines corresponding to the atoms identified for
    removal through string matching.

    Parameters:
        lines (list[str]): list of individual line strings from the CIF file.
        coordinates (list[list[float]]): list of lists containing the 3 fractional
                                        coordinates [x, y, z] of each atom due
                                        for removal.

    Returns:
        lines (list[str]): edited list of individual CIF line strings with atoms
                          at the specific coordinates removed.
        atom_count (int): total number of atoms removed from the CIF file.
    """
    atom_labels = []
    atom_count = 0

    # checking for the pymatgen formatting - different cif files
    pymatgen = lines[0]
    content = pymatgen.split(" ")
    content_fixed = [x for x in content if x]

    # setting the positions of the needed elements in the file
    if content_fixed == ["#", "generated", "using", "pymatgen\n"]:
        label_position = 1
        len_content = 6
        x = 3
        y = 4
        z = 5
    else:
        label_position = 0
        len_content = 5
        x = 2
        y = 3
        z = 4

    # going through the list of coordinates and check each line for the presense
    # of the set of three coordinates
    for coord in coordinates:
        for line in lines:
            content = line.split(" ")
            content_fixed = [x for x in content if x]
            if content_fixed[0] == "#":
                continue
            check = False

            # ignoring the lines that clearly are not coordinates
            if len(content_fixed) < len_content:
                continue

            # making the coodinares strings fit a specific pattern for comparison
            for i in range(x, z + 1):
                if content_fixed[i][0] == "1":
                    content_fixed[i] = "0" + content_fixed[i][1:]
                elif content_fixed[i][0] == "2":
                    content_fixed[i] = "1" + content_fixed[i][1:]
                elif content_fixed[i] == "-0.00000":
                    content_fixed[i] = "0.00000"
                elif content_fixed[i][0] == "-":
                    num = float(content_fixed[i])
                    if num >= -1:
                        num += 1
                    elif num >= -2:
                        num += 2
                    content_fixed[i] = format(num, ".5f")
                if len(content_fixed[i]) != 7 and content_fixed[i][0].isnumeric():
                    num = float(content_fixed[i])
                    content_fixed[i] = format(num, ".5f")
                if content_fixed[i][0].isnumeric() and "\n" in content_fixed[i]:
                    content_fixed[i] = content_fixed[i].strip()
                    num = float(content_fixed[i])
                    content_fixed[i] = format(num, ".5f")

            # comparing the coordinates to the line
            if content_fixed[x][:-1] == str(coord[0][:-1]):
                if content_fixed[y][:-1] == str(coord[1][:-1]):
                    if content_fixed[z][:-1] == str(coord[2][:-1]):
                        check = True
            # if the coordinates match the line it gets removed
            # atom label is extracted from such line
            if check is True:
                label = content_fixed[label_position]
                atom_labels.append(label)
                lines.remove(line)
                atom_count += 1

    # going through the list of atom labels to remove all the lines that
    # correspond to these atom labels
    for atom in atom_labels:
        for line in lines:
            content = line.split(" ")
            content_fixed = [x for x in content if x]
            if atom in content_fixed:
                lines.remove(line)
    return lines, atom_count

