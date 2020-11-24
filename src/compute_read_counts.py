import numpy as np
import os
import pickle as pkl

from src.coptr_contig import CoPTRContig
from src.coptr_ref import ReadFilterRef
from src.print import print_info, print_warning


def compute_read_counts_from_coverage_maps(coverage_maps):
    # instantiate classes for filtering methods

    coptr_contig = CoPTRContig(5000, 5)
    rf_ref = ReadFilterRef(5000, 0.75)
    total_passing_reads = 0
    read_counts = {}
    sample_id = None
    genome_ids = set()
    for genome_id in coverage_maps:
        cm = coverage_maps[genome_id]
        sample_id = cm.sample_id

        if cm.is_assembly and cm.passed_qc():
            binned_reads = coptr_contig.construct_coverage_matrix([cm])

            lower_bound, upper_bound = coptr_contig.compute_genomewide_bounds(binned_reads)
            count = binned_reads[np.logical_and(binned_reads >= lower_bound, binned_reads <= upper_bound)].sum()
            read_counts[cm.genome_id] = count
            total_passing_reads += count
            genome_ids.add(cm.genome_id)

        elif not cm.is_assembly:
            filtered_reads, filtered_length, qc_result = rf_ref.filter_reads(cm.read_positions, cm.length)

            if qc_result.passed_qc:
                count = filtered_reads.size
                read_counts[cm.genome_id] = count
                total_passing_reads += count

            genome_ids.add(cm.genome_id)

    rel_abun = {}
    for genome_id in read_counts:
        rel_abun[genome_id] = read_counts[genome_id] / total_passing_reads

    return sample_id, rel_abun, genome_ids




def compute_read_counts(coverage_map_folder):

    rel_abundances = {}
    genome_ids = set()
    for f in sorted(os.listdir(coverage_map_folder)):
        fname, ext = os.path.splitext(f)
        if ext != ".pkl": continue
        fpath = os.path.join(coverage_map_folder, f)

        print_info("Count", "\tprocessing {}".format(f))

        with open(fpath, "rb") as file:
            coverage_maps = pkl.load(file)

            sample_id, sample_rel_abun, sample_genome_ids = compute_read_counts_from_coverage_maps(coverage_maps)

            if sample_id is not None:
                rel_abundances[sample_id] = sample_rel_abun
                genome_ids.update(sample_genome_ids)

    return rel_abundances, genome_ids