from __future__ import annotations

import os
import csv
import shutil

from multiprocessing import current_process as cpr

from .cif_manipulation import remove_solvents_from_file

# from .log_utils import get_logger


def command_line_output(
    solvent_stats: list[str], solvent_present_flag: bool, keep_bound: bool
) -> None:
    """
    Prints information to the terminal.

    Parameters:
        solvent_stats (list[str]): list of solvent atom labels.
        solvent_present_flag (bool): True if solvent is present.
        keep_bound (bool): True if argument passed in arguments (argparse)
                           Doesn't affect bound solvent removal itself.

    Returns:
        None
    """
    if solvent_present_flag:
        if solvent_stats.get("free_solvent_flag"):
            free_solvent = solvent_stats.get("free_solvents_output")
            print(f"Identified free solvent:\n{free_solvent}")
            print("-----")
            print(f"Number of molecules: {len(free_solvent)}")
            print("-----")

        # part of output only for bound solvent
        if not keep_bound:
            if solvent_stats.get("counterions_flag"):
                counterions = solvent_stats.get("counterions_output")
                print(f"Identified counterions:\n{counterions}")
                print("-----")
                print(f"Number of molecules: {len(counterions)}")
                print("-----")
            if solvent_stats.get("bound_solvent_flag"):
                bound_solvent = solvent_stats.get("solvents_for_output")
                print(f"Identified bound solvent molecules:\n{bound_solvent}")
                print("-----")
                print(f"Number of molecules: {len(bound_solvent)}")
                print("-----")
            print("*********************")
    else:
        print("No solvent or counterions identified")
        print("*********************")


def output_csv(
    solvent_stats: list[str],
    solvent_present_flag: bool,
    total_solv_atoms: int,
    file: str,
    removed_atoms: int,
    keep_bound: bool,
) -> list[list | int | str]:
    """
    Forms a row for the solvent statistics that is appended to output csv.
    There are two types of rows: for all solvent removed and for only free solvent removed.

    Parameters:
        solvent_stats: (list[str]): list of atom labels of solvent.
        solvent_present_flag (bool): True if solvent was detected.
        total_solv_atoms (int): total number of solvent atoms.
        file (str): filename in format "*.cif".
        removed_atoms (int): number of atoms actually removed while parsing file.
        keep_bound (bool): True if argument passed in arguments (argparse)
                           Doesn't affect/change bound solvent removal itself.

    Returns:
        output_row (list[list | int | str]): formatted list of items for the csv cells.
    """
    file = os.path.basename(file)

    # forming stats for free + bound solvent export
    if keep_bound:
        # forming output row depending on presence of solvent
        if solvent_present_flag:
            # getting items from storage dict
            free_solvents_output = solvent_stats.get("free_solvents_output")
            counterions_output = solvent_stats.get("counterions_output")

            atom_count_match = True
            if total_solv_atoms != removed_atoms:
                atom_count_match = False

            # forming output row to append to the csv
            output_row = [
                file,
                "YES",
                free_solvents_output,
                len(free_solvents_output),
                counterions_output,
                len(counterions_output),
                solvent_stats.get("charge_removed"),
                total_solv_atoms,
                removed_atoms,
                atom_count_match,
                solvent_stats.get("metal_counterion_flag"),
                solvent_stats.get("huge_counterion_flag"),
            ]

        else:
            output_row = [
                file,
                "NO",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
            ]
    # forming output for only free solvent
    else:
        if solvent_present_flag:
            # getting items from storage dict
            solvents_for_output = solvent_stats.get("solvents_for_output")
            free_solvents_output = solvent_stats.get("free_solvents_output")
            counterions_output = solvent_stats.get("counterions_output")
            terminal_oxo = solvent_stats.get("oxo_mols")

            atom_count_match = True
            if total_solv_atoms != removed_atoms:
                atom_count_match = False

            # forming output row to append to the csv
            output_row = [
                file,
                "YES",
                solvents_for_output,
                len(solvents_for_output),
                free_solvents_output,
                len(free_solvents_output),
                counterions_output,
                len(counterions_output),
                terminal_oxo,
                len(terminal_oxo),
                solvent_stats.get("charge_removed"),
                total_solv_atoms,
                removed_atoms,
                atom_count_match,
                solvent_stats.get("flag_double"),
                solvent_stats.get("flag_aromatic"),
                solvent_stats.get("metal_counterion_flag"),
                solvent_stats.get("terminal_oxo_flag"),
                solvent_stats.get("entry_oxo"),
                solvent_stats.get("huge_counterion_flag"),
                solvent_stats.get("OH_removed"),
                solvent_stats.get("oxo_OH"),
            ]

        else:
            # getting the output items from the storage dictionary
            terminal_oxo_flag = solvent_stats.get("terminal_oxo_flag")
            entry_oxo = solvent_stats.get("entry_oxo")
            oxo_OH = solvent_stats.get("oxo_OH")

            # forming output row to append to the csv
            output_row = [
                file,
                "NO",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                ".",
                terminal_oxo_flag,
                entry_oxo,
                ".",
                ".",
                oxo_OH,
            ]

    return output_row


def output_cif(
    path: str, file: str, solvent_coordinates: list[list[float]], cwd: str
) -> int:
    """
    Outputs a new cif file without solvent atoms if solvent was detected.

    Parameters:
        path(str): path to the output folder.
        file (str): filename in format "*.cif".
        solvent_coordinates (list[list[float]]): list of lists containing the 3
                                        fractional coordinates [x, y, z] of each
                                        atom due for removal.
        cwd: (str): path to current working directory.

    Returns:
        removed_atoms (int): number of atoms that were removed from .cif file by the parser.
    """
    output_directory = os.path.join(path, "MOFs_removed_solvent")
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    new_cif_file = shutil.copy2(os.path.join(cwd, file), output_directory)
    with open(new_cif_file, "r+") as cif_file:
        lines = cif_file.readlines()
        cif_file.seek(0)
        cif_file.truncate()
        file_content, removed_atoms = remove_solvents_from_file(
            lines, solvent_coordinates
        )
        for line in file_content:
            cif_file.write(line)
        cif_file.close()

    return removed_atoms


def export_res(res: list[list | int | str], keep_bound: bool, output_dir) -> None:
    """
    Processes the results of the Pool to generate output csv.
    Result is a row of information/stats for one structure to be sent to the output csv.

    Parameters:
        res (list[list | int | str]): row of output information and flags for the csv.
        keep_bound (bool): True if argument passed to argparse - doesn't remove bound solvent
        output_dir (str): path to output directory

    Returns:
        None
    """
    if keep_bound:
        export_df_path = os.path.join(output_dir, "Free_solvent_removal_results.csv")
        columns = [
            "CIF",
            "Solvent",
            "Free_solvent",
            "Number_of_free_solvent_molecules",
            "Counterions",
            "Number_of_counterions",
            "Charge_removed",
            "Total_atoms",
            "Atoms_removed",
            "Atoms_match_flag",
            "Metal_counterion_flag",
            "Huge_counterion_flag",
        ]
    else:
        export_df_path = os.path.join(output_dir, "Solvent_removal_results.csv")
        columns = [
            "CIF",
            "Solvent",
            "Bound_solvent",
            "Number_of_bound molecules",
            "Free_solvent",
            "Number_of_free_solvent_molecules",
            "Counterions",
            "Number_of_counterions",
            "Terminal_oxo",
            "Number_of_terminal_oxo",
            "Charge_removed",
            "Total_atoms",
            "Atoms_removed",
            "Atoms_match_flag",
            "Flag_double",
            "Flag_aromatic",
            "Metal_counterion_flag",
            "Terminal_oxo_flag",
            "Entry_terminal_oxo",
            "Huge_counterion_flag",
            "OH_removed",
            "Oxo_OH",
        ]

    if os.path.exists(export_df_path):
        with open(export_df_path, "a", newline="") as file_obj:
            writerObj = csv.writer(file_obj)
            writerObj.writerow(res)
    else:
        with open(export_df_path, "w", newline="") as fileObj:
            writerObj = csv.writer(fileObj)
            writerObj.writerow(columns)
            writerObj.writerow(res)

