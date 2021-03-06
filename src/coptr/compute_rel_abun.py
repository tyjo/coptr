"""
compute_rel_abun.py
======================
Estimate relative abundances.
"""

"""
This file is part of CoPTR.

CoPTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CoPTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CoPTR.  If not, see <https://www.gnu.org/licenses/>.
"""

import logging
import os
import pickle as pkl

import numpy as np

from .coptr_contig import CoPTRContig
from .coptr_ref import ReadFilterRef


logger = logging.getLogger(__name__)


def compute_rel_abun_from_coverage_maps(coverage_maps, min_reads, min_cov, min_samples):
    # instantiate classes for filtering methods

    coptr_contig = CoPTRContig(min_reads, min_samples)
    rf_ref = ReadFilterRef(min_reads, min_cov)
    total_passing_reads = 0
    read_counts = {}
    genome_lengths = {}
    sample_id = None
    genome_ids = set()
    for genome_id in coverage_maps:
        cm = coverage_maps[genome_id]
        sample_id = cm.sample_id

        if cm.is_assembly and cm.passed_qc():
            binned_reads = coptr_contig.construct_coverage_matrix([cm])

            lower_bound, upper_bound = coptr_contig.compute_genomewide_bounds(
                binned_reads
            )
            count = binned_reads[
                np.logical_and(binned_reads >= lower_bound, binned_reads <= upper_bound)
            ].sum()
            read_counts[cm.genome_id] = count
            genome_lengths[cm.genome_id] = binned_reads.shape[0] * cm.compute_bin_size()
            total_passing_reads += count
            genome_ids.add(cm.genome_id)

        elif not cm.is_assembly:
            filtered_reads, filtered_length, qc_result = rf_ref.filter_reads(
                cm.read_positions, cm.length
            )

            if qc_result.passed_qc:
                count = filtered_reads.size
                read_counts[cm.genome_id] = count
                genome_lengths[cm.genome_id] = filtered_length
                total_passing_reads += count
                genome_ids.add(cm.genome_id)

    rel_abun = {}
    normalizing_constant = 0
    for genome_id in read_counts:
        unnormalized_rel_abun = (
            read_counts[genome_id] / total_passing_reads
        ) / genome_lengths[genome_id]
        rel_abun[genome_id] = unnormalized_rel_abun
        normalizing_constant += unnormalized_rel_abun

    for genome_id in read_counts:
        rel_abun[genome_id] = rel_abun[genome_id] / normalizing_constant

    return sample_id, rel_abun, genome_ids


def compute_rel_abun(coverage_map_folder, min_reads, min_cov, min_samples):

    rel_abun = {}
    genome_ids = set()
    for f in sorted(os.listdir(coverage_map_folder)):
        fname, ext = os.path.splitext(f)
        if ext != ".pkl":
            continue
        fpath = os.path.join(coverage_map_folder, f)

        logger.info("\t%s", f)

        with open(fpath, "rb") as file:
            coverage_maps = pkl.load(file)

            (
                sample_id,
                sample_rel_abun,
                sample_genome_ids,
            ) = compute_rel_abun_from_coverage_maps(
                coverage_maps, min_reads, min_cov, min_samples
            )

            if sample_id is not None:
                rel_abun[sample_id] = sample_rel_abun
                genome_ids.update(sample_genome_ids)

    return rel_abun, genome_ids
