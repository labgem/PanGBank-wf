#!/usr/bin/env python


"""Provide a command line tool to parse genome paths and taxonomy."""


import argparse
import logging
import sys
from pathlib import Path
import gzip
import scipy.sparse
from scipy.sparse import dok_matrix
from typing import Dict

from tqdm import tqdm

from scipy.cluster.hierarchy import linkage, fcluster



def remove_close_genomes(genomes_under_selection, genome_to_index, mash_dist_matrix, genomes_removed, identity_cutoff):
    """
    """
    # Get last element (which is the 'best' genome), and remove it from the list
    ref_name = genomes_under_selection.pop()

    ref_num = genome_to_index[ref_name]

    others = genomes_under_selection[::-1]

    # For each genome, compare its distance to reference genome 'ref_name'
    for gname in others:
        # Column of genome in mat_sp
        other_num = genome_to_index[gname]
        # Get distance between the 2 genomes
        if  ref_num < other_num:
            ident = mash_dist_matrix[ref_num, other_num]
        else:
            print("Should never happen as mat_sp is a triangle matrix!")
            ident = mash_dist_matrix[other_num, ref_num]

        # If distance not in the limits, remove genome from to_try and add to genomes_removed list
        if ident > identity_cutoff:
            genomes_under_selection.remove(gname)
            genomes_removed[gname] = [ref_name, identity_cutoff]

    return 0


def parse_mash_dist_result_into_matrix(genome_to_index:Dict[str, int], mash_result_file:Path, disable_bar:bool):
    """
    """

    if not mash_result_file.is_file():
        print(f"Matrix file {mash_result_file} does not exist. We cannot read it "
                     "and do the next steps. Program ending.")
        sys.exit(1)

    genome_count = len(genome_to_index)

    # Create square matrix with nbgen cols/lines. dok format is a 'Dictionary Of Keys'
    # -> writes (0, 1) value
    sparse_matrix_mash = dok_matrix((genome_count, genome_count), dtype=float)
    # Write matrix values
    proper_open = gzip.open if mash_result_file.suffix == ".gz" else open

    with tqdm(unit="k genome pair", disable=disable_bar) as progress:
        with proper_open(mash_result_file, "rt") as matf:
            for i, line in enumerate(matf):
                path1, path2, dist = line.split()[:3]
                num1 = genome_to_index[path1]
                num2 = genome_to_index[path2]

                # only in lower triangle (no duplicate)
                if  num1 == num2:
                    pass
                elif num1 < num2:
                    sparse_matrix_mash[num1, num2] = 1 - float(dist)
                else:
                    sparse_matrix_mash[num2, num1] = 1 - float(dist)

                if i % 10000 == 0:
                    progress.update(10)

    return sparse_matrix_mash

def select_genomes_like_panacota(mash_dist_matrix, sorted_genomes, genome_to_index, identity_cutoff, disable_bar):
    genomes_under_selection = sorted_genomes[::-1]
    # Put list of genomes removed by mash comparison, and why
    # (out of limits distance with which genome)
    genomes_removed = {}  # {genome: [compared_with, dist]}
    genome_count = len(genomes_under_selection)

    last_genome_process_count = 0

    with tqdm(total=genome_count, unit="genome", disable=disable_bar) as progress:

        while len(genomes_under_selection) > 1:

            remove_close_genomes(genomes_under_selection, genome_to_index, mash_dist_matrix, genomes_removed, identity_cutoff)

            current_genome_process_count = genome_count - len(genomes_under_selection)
            genome_process_in_loop = current_genome_process_count - last_genome_process_count
            last_genome_process_count = current_genome_process_count
            progress.update(genome_process_in_loop)



    logging.info("Final number of genomes in dataset: {}".format(genome_count - len(genomes_removed)))

    selected_genomes = [genome for genome in sorted_genomes if genome not in genomes_removed]


    return selected_genomes

def select_genome_with_hierarchical_clustering(mash_dist_matrix, num_clusters, index_to_genome, method):

    Z = linkage(mash_dist_matrix, method=method)

    # Use fcluster to cut the dendrogram and get cluster labels
    cluster_labels = fcluster(Z, num_clusters, criterion='maxclust')

    known_clusters = set()
    cluster_to_genomes = {}
    selected_genomes = []
    for i, cluster in enumerate(cluster_labels):
        if cluster not in known_clusters:
            # print(f"genome {i} in cluster {cluster}")
            selected_genomes.append(index_to_genome[i])

            known_clusters.add(cluster)
            cluster_to_genomes[cluster] = [i]
        else:
            cluster_to_genomes[cluster].append(i)


    # for cluster, genomes in cluster_to_genomes.items():
    #     print(f"cluster={cluster} (size={len(genomes)}): {' '.join(map(str, genomes[:20]))}...")

    print(len(selected_genomes), "Cluster")

    return selected_genomes


def maxmin_clustering_with_sorted_quality(dist_matrix, n):
    """
    Selects a subset of n genomes using the MaxMin algorithm, starting
    with the highest quality genome (index 0).

    Parameters:
    - dist_matrix: np.array of distances (M x M)
    - n: Number of genomes to select

    Returns:
    - selected_indices: Indices of the selected genomes
    """
    num_genomes = dist_matrix.shape[0]

    # Start with the first genome (best quality, index 0)
    selected_indices = [0]

    # Iterate until n genomes are selected
    while len(selected_indices) < n:
        max_dist = -1
        next_idx = -1

        # Iterate over the remaining genomes (already sorted by quality)
        for i in range(num_genomes):
            if i not in selected_indices:
                # Calculate the minimum distance between this genome
                # and the already selected genomes
                min_dist_to_selected = min([dist_matrix[i, j] for j in selected_indices])

                # Select the genome that maximizes this minimum distance
                if min_dist_to_selected > max_dist:
                    max_dist = min_dist_to_selected
                    next_idx = i

        # Add the most distant genome to the subset
        selected_indices.append(next_idx)

    return selected_indices


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse genome paths and taxonomy.",
        epilog="Example: python parse_genomes_and_taxonomy.py --genomes <genomes_file> --taxonomy <taxonomy_file>",
    )
    parser.add_argument(
        "--sorted_genomes_file",
        type=Path,
        required=True,
        help="Path to a file listing genome path file sorted from best to worst "
    )
    parser.add_argument(
        "--mash_dist_result",
        type=Path,
        required=True,
        help="Result file of mash dist",
    )

    parser.add_argument(
        "-d", "--min_dist",
        help="Discard genome(s) closer than a Mash distance. ",
        default=0.005,
        type=float,
    )

    parser.add_argument(
        "--disable_bar",
        type=bool,
        default=False,
        help="Disable progress bar",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Directory where diltered genomes are stored.",
        default="genome_derep_out",
        type=Path,
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="Desired log level (default: WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )

    return parser.parse_args(argv)


def write_selected_genomes_ids(selected_genomes, outfile):


    with open(outfile, "w") as fl:
        fl.write("\n".join(selected_genomes))


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.sorted_genomes_file.is_file():
        logging.error(f"Sorted genome file {args.sorted_genomes_file} was not found!")
        sys.exit(2)

    if not args.mash_dist_result.is_file():
        logging.error(f"mash_dist_result file {args.mash_dist_result} was not found!")
        sys.exit(2)

    if not args.output.is_dir():
        logging.debug(f"Create output directory {args.output.absolute().as_posix()}")
        Path.mkdir(args.output, exist_ok=True)


    sorted_genomes_file = args.sorted_genomes_file
    with open(sorted_genomes_file) as fl:
        sorted_genomes = [line.rstrip() for line in fl]


    genome_to_index = {genome: index for index, genome in enumerate(sorted_genomes)}

    sparse_matrix_file:Path = args.output / "sparse_matrix_mash_dist.npz"

    if sparse_matrix_file.is_file():
        logging.info(f"Loading matrix contained in {sparse_matrix_file}")
        # convert matrix returned by load_npz (coo format, as saved) to dok format
        sparse_similarity_matrix = scipy.sparse.load_npz(sparse_matrix_file).todok()

    else:
        sparse_similarity_matrix = parse_mash_dist_result_into_matrix(genome_to_index, args.mash_dist_result, disable_bar=args.disable_bar)

        logging.info("Saving matrix to npz file to be loaded quicker if needed later")
        # Convert dok_matrix to coo format, as dok format is not allowed by save_npz
        coo_mat = sparse_similarity_matrix.tocoo()
        scipy.sparse.save_npz(sparse_matrix_file, coo_mat)

    for distance_cutoff in [0, 0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007]:
        identity_cutoff = 1 - distance_cutoff

        # Panacota selection
        # start_time = time.time()
        selected_genomes_panacota = select_genomes_like_panacota(sparse_similarity_matrix, sorted_genomes, genome_to_index, identity_cutoff, args.disable_bar)
        # panacota_time = time.time() - start_time

        selected_genome_dir = args.output / str(distance_cutoff)
        selected_genome_outfile = selected_genome_dir / "selected_genomes.list"
        Path.mkdir(selected_genome_dir, exist_ok=True)

        write_selected_genomes_ids(selected_genomes_panacota, outfile=selected_genome_outfile)



if __name__ == "__main__":
    main()
