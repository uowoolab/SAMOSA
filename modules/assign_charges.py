from __future__ import annotations

# handle differing python versions (< 3.8, & > 3.9) of CSD python API
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

import mendeleev
from ccdc.molecule import Atom, Molecule
from ccdc.descriptors import MolecularDescriptors

from collections import defaultdict
from multiprocessing import current_process as cpr

from .log_utils import get_logger


def get_unique_sites(mole: Molecule, asymmole: Molecule) -> list[Atom]:
    """
    Get the unique atoms in a MOF structure that belong to the asymmetric unit.

    Parameters:
        mole (ccdc.molecule.Molecule): The MOF structure to get unique atoms from.
        asymmole (ccdc.molecule.Molecule): The asymmetric unit of the MOF structure.

    Returns:
        uniquesites (list[ccdc.molecule.Atom]): list of unique atoms in the MOF structure
                                                that belong to the asymmetric unit.
    """
    # blank list for unique sites
    uniquesites = []
    labels = []
    asymmcoords = []
    molecoords = []
    duplicates = []
    for atom in asymmole.atoms:
        asymmcoords.append(
            atom.coordinates
        )  # writing coordinates of atoms in asymmetric unit
    for atom in mole.atoms:
        if (
            atom.coordinates in asymmcoords
        ):  # if they are coordinates of atom in asymmetric unit
            if (
                atom.coordinates not in molecoords
            ):  # and they are not in list of coordinates in whole MOF
                if atom.label not in labels:  # and we didn't have this atom before
                    uniquesites.append(atom)
                    molecoords.append(atom.coordinates)
                    labels.append(atom.label)
                else:
                    duplicates.append(atom)
            else:
                duplicates.append(atom)
    if len(duplicates) >= 1:
        for datom in duplicates:
            for atom in uniquesites:
                if any(
                    [
                        (datom.coordinates == atom.coordinates),
                        (datom.label == atom.label),
                    ]
                ):
                    if datom.atomic_symbol == atom.atomic_symbol:
                        if len(datom.neighbours) > len(atom.neighbours):
                            uniquesites.remove(atom)
                            uniquesites.append(datom)
                    elif datom.label not in labels:
                        uniquesites.append(datom)
    return uniquesites


def get_metal_sites(sites: list[Atom]) -> list[Atom]:
    """
    Get the metal sites in a MOF structure that belong to the asymmetric unit.

    Parameters:
        sites (list[ccdc.molecule.Atom]):  list of unique atoms in the MOF structure
                                           that belong to the asymmetric unit.

    Returns:
        metalsites (list[ccdc.molecule.Atom]): list of metal sites in the MOF structure
                                               that belong to the asymmetric unit.
    """
    metalsites = []
    for site in sites:
        if site.is_metal is True:
            metalsites.append(site)
    return metalsites


def get_binding_sites(
    metalsites: list[Atom], uniquesites: list[Atom]
) -> tuple[list[Atom], dict[Atom, list[Atom]]]:
    """
    Get the binding sites in a MOF structure, given the list of unique metal atoms
    and all unique atoms in the MOF.

    Parameters:
        metalsites (list[ccdc.molecule.Atom]): list of unique metal atoms in the MOF.
        uniquesites (list[ccdc.molecule.Atom]): list of unique atoms in the MOF.

    Returns:
        binding_sites (list[ccdc.molecule.Atom]): list of binding sites connecting metal
                                                  atoms and ligands.
        binding_pairs (dict[ccdc.molecule.Atom, list[ccdc.molecule.Atom]]): dictionary
                                with each ligand binding site Atom object as keys and
                                a list of the associated bound metal atoms as values.
    """
    binding_sites = []
    binding_pairs = defaultdict(list)
    for metal in metalsites:
        for ligand in metal.neighbours:
            for site in uniquesites:
                if ligand.label == site.label:
                    binding_sites.append(site)
                    binding_pairs[site].append(metal)
    return binding_sites, binding_pairs


def ringVBOs(mole: Molecule) -> dict[int, int]:
    """
    Calculates the VBO (valence bond order) for each atom in the MOF.

    Parameters:
        mole (ccdc.molecule.Molecule): molecule object of MOF structure
    Returns:
        ringVBO (dict[int, int]): dictionary with each atom's index in mole.atoms
                                 as keys and VBO (valence bond order) as values
    """
    ringVBO = {}
    unassigned = mole.atoms
    ringcopy = mole.copy()
    oncycle_atoms = []
    offcycle_atoms = []
    oncycle_labels = []
    offcycle_labels = []
    # remove all the metals, this
    # prevents metal-containing rings (i.e. pores)
    # from interfering
    for atom in ringcopy.atoms:
        if atom.is_metal:
            ringcopy.remove_atom(atom)
    # collect all the cyclic atoms
    for atom in ringcopy.atoms:
        if atom.is_cyclic:
            if atom not in oncycle_atoms:
                oncycle_atoms.append(atom)
                oncycle_labels.append(atom.label)
    # we also need everything that the cyclic atoms are bound to
    for atom in oncycle_atoms:
        for neighbour in atom.neighbours:
            if neighbour not in oncycle_atoms:
                if neighbour not in offcycle_atoms:
                    offcycle_atoms.append(neighbour)
                    offcycle_labels.append(neighbour.label)
    cyclicsystem = oncycle_atoms + offcycle_atoms
    # remove every atom that isn't part of or directly bound to a cycle
    for atom in ringcopy.atoms:
        if atom not in cyclicsystem:
            ringcopy.remove_atom(atom)
    # find all non-cyclic bonds
    # bonds between cycles, break and cap with H
    for bond in ringcopy.bonds:
        if not bond.is_cyclic:
            # bonds between cycles
            if all((member.label in oncycle_labels for member in bond.atoms)):
                member1 = bond.atoms[0]
                member2 = bond.atoms[1]
                Hcap1 = Atom("H", coordinates=member1.coordinates)
                Hcap2 = Atom("H", coordinates=member2.coordinates)
                Hcap1_id = ringcopy.add_atom(Hcap1)
                Hcap2_id = ringcopy.add_atom(Hcap2)
                ringcopy.add_bond(bond.bond_type, Hcap1_id, member2)
                ringcopy.add_bond(bond.bond_type, Hcap2_id, member1)
                ringcopy.remove_bond(bond)
    # cap off-cycle atoms
    for offatom in offcycle_atoms:
        # get the VBO for each off-cycle atom
        # (VBO with respect to cyclic atoms)
        offVBO = 0
        # quick check for delocalised systems in the ring
        # if there are any, get the delocalised bond orders
        if any(bond.bond_type == "Delocalised" for bond in offatom.bonds):
            offdVBO = delocalisedLBO(offcycle_atoms)
        for bond in offatom.bonds:
            # Each bond contributes to Ligand Bond Order according to its type
            if bond.bond_type == "Single":
                offVBO += 1
            elif bond.bond_type == "Double":
                offVBO += 2
            elif bond.bond_type == "Triple":
                offVBO += 3
            elif bond.bond_type == "Quadruple":
                offVBO += 4
            elif bond.bond_type == "Delocalised":
                offVBO += offdVBO[offatom]
            elif bond.bond_type == "Aromatic":
                offVBO += 0
                get_logger(cpr().name).warning("impossible aromatic bond detected")
        # cap with appropriate element for VBO
        if offVBO == 1:
            offatom.atomic_symbol = "H"
        elif offVBO == 2:
            offatom.atomic_symbol = "O"
        elif offVBO == 3:
            offatom.atomic_symbol = "N"
        elif offVBO == 4:
            offatom.atomic_symbol = "C"
        elif offVBO == 5:
            offatom.atomic_symbol = "P"
        elif offVBO == 6:
            offatom.atomic_symbol = "S"
        elif offVBO > 6:
            get_logger(cpr().name).warning(
                "issue detected in valence bond order calculations (capping)"
            )
    # for each cyclic system, reassign bonds, kekulize, and get VBO
    # the bond and atom pruning we did above ensures that fused cycles
    # will be treated as a single system
    # while non-fused cycles that are connected via bonding are treated
    # as seperate systems
    for cyclesys in ringcopy.components:
        # reassign bonds and kekulize
        cyclesys.assign_bond_types()
        cyclesys.kekulize()
        # quick check for delocalized systems in the ring
        # if there are any, get the delocalised bond orders
        if any(bond.bond_type == "Delocalised" for bond in cyclesys.bonds):
            rdVBO = delocalisedLBO(cyclesys)
        # assign VBO for each on-cycle atom
        for ratom in cyclesys.atoms:
            rVBO = 0
            if ratom.label in oncycle_labels:
                for rbond in ratom.bonds:
                    # Each bond contributes to Ligand Bond Order according to its type
                    if rbond.bond_type == "Single":
                        rVBO += 1
                    elif rbond.bond_type == "Double":
                        rVBO += 2
                    elif rbond.bond_type == "Triple":
                        rVBO += 3
                    elif rbond.bond_type == "Quadruple":
                        rVBO += 4
                    elif rbond.bond_type == "Delocalised":
                        rVBO += rdVBO[ratom]
                    elif rbond.bond_type == "Aromatic":
                        rVBO += 0
                        get_logger(cpr().name).warning(
                            "impossible aromatic bond detected"
                        )
                # the VBOs are currently associated to atom objects
                # in molecule objects that we have modified
                # we need these to be associated to atom objects in
                # the parent (unmodified) molecule object
                for matom in unassigned:
                    if matom.label == ratom.label:
                        ringVBO[matom] = rVBO
                        unassigned.remove(matom)
    return ringVBO


def assign_VBS(atom: Atom, rVBO: dict[int, int], dVBO: dict[int, float]) -> int:
    """
    Assigns a Valence-Bond-Sum (VBS) to an atom.

    Parameters:
        atom (ccdc.molecule.Atom): atom object
        rVBO (dict[int, int]): dictionary with each atom's index in mole.atoms
                              as keys and VBO (valence bond order) as values
        dVBO (dict[int, float]): dictionary with delocalized bond-possessing atom's
                                index in mole.atoms as keys and their corresponding
                                (delocalized-only) VBS.

    Returns:
        VBO (int): valence bond sum value.
    """
    VBO = 0
    if atom.is_metal:
        return 0
    if atom in rVBO:
        VBO = rVBO[atom]
    else:
        for bond in atom.bonds:
            if any(batom.is_metal for batom in bond.atoms):
                VBO += 0
            # Each bond contributes to Ligand Bond Order according to its type
            elif bond.bond_type == "Single":
                VBO += 1
            elif bond.bond_type == "Double":
                VBO += 2
            elif bond.bond_type == "Triple":
                VBO += 3
            elif bond.bond_type == "Quadruple":
                VBO += 4
            elif bond.bond_type == "Delocalised":
                VBO += dVBO[atom]
            elif bond.bond_type == "Aromatic":
                VBO += rVBO[atom]
    return VBO


def delocalisedLBO(molecule: Molecule) -> dict[int, float]:
    """
    Takes a ccdc.molecule.Molecule as input and writes a dictionary of all atoms
    in the molecule with delocalized bonds and their (delocalized-only) VBS.

    Parameters:
        molecule (ccdc.molecule.Molecule):  MOF molecule object.

    Returns:
        delocal_dict (dict[int, float]): dictionary with delocalized bond-possessing
                                        atom's index in mole.atoms as keys and their
                                        corresponding (delocalized-only) VBS.
    """

    def TerminusCounter(atomlist: list[Atom]) -> int:
        """
        Counts the number of termini in the input delocalized bond system.

        Parameters:
            atomlist (list[ccdc.molecule.Atom]): list of atoms in delocalised system.

        Returns:
            NTerminus (int): number of termini in delocalized bond system.
        """
        NTerminus = 0
        for member in atomlist:
            connectivity = 0
            for bond in member.bonds:
                if bond.bond_type == "Delocalised":
                    connectivity += 1
            if connectivity == 1:
                NTerminus += 1
        return NTerminus

    def delocal_crawl(atomlist: list[Atom]) -> list[Atom]:
        """
        Recursively searches for atoms in delocalised bond systems starting from
        an input list containing at least one delocalised bonding atom.

        Parameters:
            atomlist (list[ccdc.molecule.Atom)]: list of atoms in delocalised system.

        Returns:
            atomlist (list[ccdc.molecule.Atom]): modified list of atoms in
                                               delocalised system.
        """
        for delocatom in atomlist:
            for bond in delocatom.bonds:
                if bond.bond_type == "Delocalised":
                    for member in bond.atoms:
                        if member not in atomlist:
                            atomlist.append(member)
                            return delocal_crawl(atomlist)
        return atomlist

    delocal_dict = {}
    # Handle cases where ccdc.Molecule() or atom list is input
    molecule = molecule if isinstance(molecule, list) else molecule.atoms
    for atom in molecule:
        if all(
            [
                (any(bond.bond_type == "Delocalised" for bond in atom.bonds)),
                (atom not in delocal_dict),
            ]
        ):
            delocal_dict[atom] = []
            delocal_system = delocal_crawl([atom])
            NTerminus = TerminusCounter(delocal_system)
            for datom in delocal_system:
                connectivity = 0
                delocLBO = 0
                for neighbour in datom.neighbours:
                    if neighbour in delocal_system:
                        connectivity += 1
                if connectivity == 1:
                    # terminus
                    delocLBO = (NTerminus + 1) / NTerminus
                if connectivity > 1:
                    # node
                    delocLBO = (connectivity + 1) / connectivity
                delocal_dict[datom] = delocLBO
    return delocal_dict


def get_CN(atom: Atom) -> int:
    """
    Returns the coordination number of an atom.

    Parameters:
        atom (ccdc.molecule.Atom): atom object.

    Returns:
        coord_number (int): atom's coordination number.
    """
    coord_number = 0
    for neighbour in atom.neighbours:
        if not neighbour.is_metal:
            coord_number += 1
    return coord_number


def valence_e(atom: Atom) -> int:
    """
    Returns number of valence electrons of an atom.

    Parameters:
        atom (ccdc.molecule.Atom): atom object.

    Returns:
        valence (int): atom's valence electron count.
    """
    elmnt = mendeleev.element(atom.atomic_symbol)
    if elmnt.block == "s":
        valence = elmnt.group_id
    elif elmnt.block == "p":
        valence = elmnt.group_id - 10
    elif elmnt.block == "d":
        valence = elmnt.group_id
    elif elmnt.block == "f":
        if elmnt.atomic_number in range(56, 72):
            valence = elmnt.atomic_number - 57 + 3
        elif elmnt.atomic_number in range(88, 104):
            valence = elmnt.atomic_number - 89 + 3
        else:
            get_logger(cpr().name).error("unexpected f block element")
            raise ValueError("valence_e() >> Unexpected f block element", elmnt)
    elif elmnt.group_id == 18:
        valence = 8 if elmnt.symbol != "He" else 2
    else:
        get_logger(cpr().name).error(
            "unexpected element in valence electron calculations"
        )
        raise ValueError("valence_e() >> Unexpected valence electrons", elmnt)
    return valence


def carbocation_check(atom: Atom) -> Literal["tetrahedral", "trigonal"]:
    """
    Check carbocation/carbanion geometry.

    Parameters:
        atom (ccdc.molecule.Atom): atom object.

    Returns:
        Literal["tetrahedral", "trigonal"]: geometry at input atom.
    """
    abc = []
    # get atom neighbours
    for neighbours in atom.neighbours:
        if not neighbours.is_metal:
            abc.append(neighbours)
    # get all three relevant bond angles
    angle1 = MolecularDescriptors.atom_angle(abc[0], atom, abc[1])
    angle2 = MolecularDescriptors.atom_angle(abc[0], atom, abc[2])
    angle3 = MolecularDescriptors.atom_angle(abc[1], atom, abc[2])
    # average the angels
    AVGangle = abs(angle1 + angle2 + angle3) / 3
    # take the difference between the averaged bond angles and
    # ideal trigonal planar/tetrahedral bond angles
    tet = abs(AVGangle - 109.5)
    trig = abs(AVGangle - 120)
    if tet < trig:
        return "tetrahedral"
    if trig < tet:
        return "trigonal"


def carbene_type(atom: Atom) -> Literal["singlet", "triplet"]:
    """
    Distinguishes between singlet and triplet carbenes.

    Parameters:
        atom (ccdc.molecule.Atom): atom object suspected of being carbene
                                  (2-coordinate carbon II)

    Returns:
        Literal["singlet", "triplet"]: carbene type at input atom.
    """
    # get alpha-atoms
    alpha = atom.neighbours
    alpha_type = []
    # get element symbols for alpha atoms
    for a in alpha:
        if not a.is_metal:
            alpha_type.append(a.atomic_symbol)
    # if any alpha atom is a heteroatom, return "singlet"
    # these are Fischer carbenes
    for a in alpha_type:
        if not any([(a == "C"), (a == "H")]):
            return "singlet"
    # if the carbene C is in a heterocycle,
    # return "singlet"
    # there are Arduengo carbenes (NHCs, CAACs)
    if atom.is_cyclic is True:
        for ring in atom.rings:
            for species in ring.atoms:
                if not species.atomic_symbol == "C":
                    return "singlet"
    # for all other carbenes, return "triplet"
    # these are Schrock carbenes
    return "triplet"


def iVBS_Oxidation_Contrib(
    unique_atoms: list[Atom], rVBO: dict[int, int], dVBO: dict[int, float]
) -> dict[Atom, float]:
    """
    Determines the oxidation state contribution for all unique atoms.

    Parameters:
        unique_atoms (list[ccdc.molecule.Atom]): unique atoms belonging to the
                                                asymmetric unit.
        rVBO (dict[int, int]): dictionary with each atom's index in mole.atoms
                              as keys and VBO (valence bond order) as values
        dVBO (dict[int, float]): dictionary with delocalized bond-possessing atom's
                                index in mole.atoms as keys and their corresponding
                                (delocalized-only) VBS.

    Returns:
        oxi_contrib (dict[ccdc.molecule.Atom, float)]: dictionary with Atom object as keys
                                                      and their oxidation state contribution
                                                      as values
    """
    VBS = 0
    CN = 0
    valence = 0
    oxi_contrib = {}
    # for each unique atom
    for atom in unique_atoms:
        # assign valence-bond-sum
        VBS = assign_VBS(atom, rVBO, dVBO)
        # determine coordination number
        CN = get_CN(atom)
        #  determine number of valence electrons
        valence = valence_e(atom)
        # get number of unpaired electrons in the free element
        unpaired_e = 4 - abs(4 - valence)

        #  metals do not contribute:
        if atom.is_metal:
            oxi_contrib[atom] = 0
        # Normal valences:
        elif VBS <= (unpaired_e):
            oxi_contrib[atom] = unpaired_e - VBS
        # Expanded (2e) valences:
        elif (VBS > unpaired_e) and (VBS < valence):
            diff = VBS - unpaired_e
            if diff <= 2:
                UPE = valence - unpaired_e - 2
            elif diff <= 4:
                UPE = valence - unpaired_e - 4
            elif diff <= 6:
                UPE = valence - unpaired_e - 6
            elif diff <= 8:
                UPE = valence - unpaired_e - 8
            oxi_contrib[atom] = VBS + UPE - valence
        elif VBS >= (valence):
            oxi_contrib[atom] = VBS - valence

        # need to check for 3-coordinate carbocations,
        # 3-coordinate carbanions, carbenes, and heavier
        # homologues (these are not immediately detectable)
        if any(
            [
                (atom.atomic_symbol == "C"),
                (atom.atomic_symbol == "Si"),
                (atom.atomic_symbol == "Ge"),
                (atom.atomic_symbol == "Pb"),
            ]
        ):
            if atom not in rVBO:
                # 3 coordinate and VBS 3 could be
                # carbanion or carbocation
                if VBS == 3 and CN == 3:
                    geom = carbocation_check(atom)
                    if geom == "trigonal":
                        oxi_contrib[atom] = -1
                    if geom == "tetrahedral":
                        oxi_contrib[atom] = 1
            # VBS 2 and 2 coordinate is carbene,
            # but singlet or triplet?
            if VBS == 2 and CN == 2:
                carbene = carbene_type(atom)
                if carbene == "singlet":
                    oxi_contrib[atom] = 2
                if carbene == "triplet":
                    oxi_contrib[atom] = 0

        # Nitro groups frequently have both N-O bonds assigned
        # as double bonds, giving incorrect VBS of 5
        # and oxidation contribution of -2
        # this block catches this and applies a fix
        if all(
            [
                (atom.atomic_symbol == "N"),
                (VBS == 5 and CN == 3),
            ]
        ):
            N_sphere1 = atom.neighbours
            O_count = 0
            for neighbour in N_sphere1:
                if neighbour.atomic_symbol == "O":
                    O_count += 1
            geom = carbocation_check(atom)
            if O_count == 2 and geom == "trigonal":
                oxi_contrib[atom] = 0
    return oxi_contrib


def redundantAON(AON: dict[Atom, float], molecule: Molecule) -> dict[Atom, float]:
    """
    Takes the oxidation contributions of atoms in unique sites and assigns them
    to the atoms in redundant sites based on equality of atom labels.

    Parameters:
        AON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                                              and their oxidation state contribution
                                              as values for unique Atom objects.
        molecule (ccdc.molecule.Molecule): Molecule object.

    Returns:
        redAON (dict[ccdc.molecule.Atom, float]): dictionary with Atom object as keys
                                                 and their oxidation state contribution
                                                 as values for all (including redundant)
                                                 Atom objects.
    """
    redAON = {}
    for rsite1 in molecule.atoms:
        for usite1 in AON:
            redAON[usite1] = AON[usite1]
            if rsite1.label == usite1.label:
                redAON[rsite1] = AON[usite1]
    return redAON

