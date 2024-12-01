#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import logging
import time

# from itertools import repeat
from multiprocessing import Manager, Queue, Pool, current_process

from modules.assign_charges import *
from modules.cif_manipulation import *
from modules.outputs import *
from modules.solvent_analysis import *

from modules.log_utils import main_logging_setup, child_logging_setup


code_desc = """
Structural Activation via Metal Oxidation State Analysis.

A solvent removal protocol generating activated crystal structures from experimental crystallographic information.
Removes free solvent, counterion, and/or bound solvent molecules keeping into account the charge of removed fragments.

Designed for crystal structures in the CIF file format.
Primary intent of the method is to clean Metal-Organic Frameworks (MOFs) for computation, but it can be applied to any 3D periodic crystal structure.
"""
parser = argparse.ArgumentParser(
    description=code_desc, formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    "--files_path", default=os.getcwd(), help="specify the folder where MOFs are stored"
)
parser.add_argument(
    "--export_path", default=os.getcwd(), help="specify the output folder"
)
parser.add_argument(
    "--n_processes", default=4, help="specify number of parallel processes"
)
parser.add_argument(
    "--removable_denticity",
    type=int,
    default=1,
    help="maximum value of ligand denticity that will still be considered for removal (default: 1)",
)
parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="if set, turns on command line outputs and optional log messages",
)
parser.add_argument(
    "--keep_bound", action="store_true", help="if set, bound solvent is not removed"
)
parser.add_argument(
    "--keep_oxo", action="store_true", help="if set, terminal oxygens are not removed"
)
parser.add_argument(
    "--logging",
    default="INFO",
    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    help="Set the logging level (default: INFO).",
)
args = parser.parse_args()

f_path = args.files_path
output_dir = args.export_path
keep_bound = args.keep_bound
keep_oxo = args.keep_oxo
log_type = getattr(logging, args.logging, logging.DEBUG)


def worker(file: str, queue: Queue) -> dict[str, list | int | str] | None:
    """
    Handles the overall SAMOSA workflow including reading in the structure,
    identify and removing solvent molecules, and outputting the final structure
    file (CIF) and summary (csv). Separate function call for each structure/process
    in the multiprocessing Pool.

    Parameters:
        file (str): filename in format "*.cif".
        queue (multiprocessing.Queue): shared queue for child processes logging.

     Returns:
        output_row (dict[str, list | int | str]): collects statistics for
                                                  output csv.
    """
    worker_logger = child_logging_setup(queue, current_process().name, log_type)
    try:
        start_time_MOF = time.time()
        worker_logger.info("analyzing %s ..." % os.path.basename(file))
        input_cif = file

        worker_logger.debug("Reading structure file ...")
        # read in the .cif, extract the underlying molecule,
        # identify the unique sites, metal sites, binding sites
        cif = readentry(input_cif)
        mol = cif.molecule
        asymmol = cif.asymmetric_unit_molecule

        worker_logger.debug("Analyzing sites ...")
        uniquesites = get_unique_sites(mol, asymmol)
        metalsites = get_metal_sites(uniquesites)
        binding_sites, binding_pairs = get_binding_sites(
            metalsites, uniquesites
        )  # outputs list of atoms bound to metals

        worker_logger.debug("Calculating oxidation state contributions ...")
        # Now get the localized oxidation state contribution of each atom
        # need delocalized bond contributions
        dVBO = delocalisedLBO(mol)
        # then need ring bond contributions - NEW TO dev_194
        rVBO = ringVBOs(mol)
        # finally combine delocal/aromatic bond contributions with localized bonding
        AON = iVBS_Oxidation_Contrib(uniquesites, rVBO, dVBO)
        # Previous only assigns an oxidation contribution to unique images of atoms,
        # also need to assign these values to redundant sites:
        rAON = redundantAON(AON, mol)

        # the atoms in rAON calculated and atoms in molecule_no_metals are different
        # redoing the rAON dictionary to have atom labels instead of atoms
        rAON_atomlabels = get_rAON_atomlabels(rAON)

        worker_logger.debug("Entering free solvent removal ...")
        # creates a copy of the initial molecule
        # cleans the molecule from free solvents
        (
            molecule_work,
            free_solvents,
            counterions,
            solvent_stats_batch1,
        ) = remove_free_solvent(mol, rAON_atomlabels)

        # if the user wants to keep bound solvent skipping further steps
        if not keep_bound:
            worker_logger.debug("Entering bound solvent removal ...")
            # identifying the oxo molecules
            oxo_mols, solvent_stats_batch2 = get_oxo(uniquesites, file, keep_oxo)

            # removes metals from the molecule
            molecule_no_metals = remove_metals(molecule_work)

            # creates the initial list of possible solvents
            solvent_mols, solvent_stats_batch3 = define_solvents(
                molecule_no_metals, rAON_atomlabels
            )

            # check the possible solvents
            solvent_mols_checked, solvent_stats_batch4 = check_solvent(
                solvent_mols, binding_pairs, uniquesites, args.removable_denticity
            )

            # combine all the lists of items to remove
            solvents_to_remove = get_solvents_to_remove(
                solvent_mols_checked,
                free_solvents,
                counterions,
                oxo_mols,
            )

            # putting together dictionary of otput stats
            output_solvent_stats = {
                **solvent_stats_batch1,
                **solvent_stats_batch2,
                **solvent_stats_batch3,
                **solvent_stats_batch4,
            }
        else:
            solvents_to_remove = free_solvents + counterions
            output_solvent_stats = solvent_stats_batch1

        # checking if there is any solvent
        solvent_present_flag = False
        if len(solvents_to_remove) > 0:
            solvent_present_flag = True

        worker_logger.debug("Generating output cif ...")
        # getting the coordinates of the atoms for removal
        solvent_coordinates = get_coordinates(mol, solvents_to_remove)

        total_solv_atoms = len(solvents_to_remove)

        # removing solvent from cif file if solvent is present
        if solvent_present_flag:
            removed_atoms = output_cif(output_dir, file, solvent_coordinates, f_path)
        else:
            removed_atoms = 0

        worker_logger.debug("Generating output csv ...")
        # output a csv with solvent removal stats
        output_row = output_csv(
            output_solvent_stats,
            solvent_present_flag,
            total_solv_atoms,
            file,
            removed_atoms,
            keep_bound,
        )

        if args.verbose:
            command_line_output(output_solvent_stats, solvent_present_flag, keep_bound)

        worker_logger.info(
            "%s / summary > %s"
            % (os.path.basename(file), ",".join([str(j) for j in output_row]))
        )
        worker_logger.info(
            "%s finished in %s seconds"
            % (os.path.basename(file), (time.time() - start_time_MOF))
        )
        return output_row

    except Exception as e:
        worker_logger.error("%s failed" % file)


def main():
    start_time = time.time()
    # init manager to handle process logging queue
    with Manager() as manager:
        log_q = manager.Queue()
        root_logger, listener = main_logging_setup(log_q, log_type)
        # search for structure files
        files = [
            os.path.join(f_path, file)
            for file in os.listdir(f_path)
            if file.endswith(".cif")
        ]
        root_logger.info("Selected logging level: logging.%s" % (args.logging))
        root_logger.info("%s cif files detected in %s" % (len(files), f_path))
        # init multiprocessing pool
        with Pool(processes=int(args.n_processes)) as pool:
            # f_args = [(f, log_q) for f in files]
            for res in pool.starmap(worker, [(f, log_q) for f in files]):
                if res is not None:
                    export_res(res, keep_bound, output_dir)
        # cleanup
        root_logger.info(
            "SAMOSA finished running in %s seconds" % (time.time() - start_time)
        )
        # end queue and listener
        log_q.put_nowait(None)
        listener.stop()


if __name__ == "__main__":
    main()

