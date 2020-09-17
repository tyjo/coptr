"""
bam_processor.py
======================
Extract coordinates along each reference sequence and group reference 
sequences by genome.
"""

import math
import numpy as np
import os.path
import pysam
import re

from src.print import print_info

class BamProcessor:
    """Extract coordinates from a bam file.

    Parameters
    ----------
        regex : str
            A regular expression that matches a reference genome id from
            a reference sequence name.

    Notes
    -----
        For assemblies, the reference database will contain several contigs
        per assembly, each with a reference sequence name. We need to map
        a reference sequence name to the assembly is came from. This is accomplished
        with a regular expression. The expression matches a reference genome id.

        For example, consider an assembly with contigs named

        | ``ERS235517|65|k99_317``
        | ``ERS235517|65|k99_655``
        | ``ERS235517|65|k99_1708``

        ERS235517|65 gives the reference genome, while the k99_* identify
        each contig. We want a regular expression that matches ERS235517|65.
        The regular expression

        ``\w+\|\\d+``

        matches 1 or more letter, number, or underscore (``\w+``), followed by
        the ``|`` character, followed by 1 or more number (``\d+``). This extracts the
        right reference genome id.
    """

    def __init__(self, regex="[^\|]+"):
        self.regex = regex
        self.min_avg_qual = 30

        # bam files without an index will generate a warning on
        # opening. since we don't need an index, setting the
        # verbosity will silence this message
        pysam.set_verbosity(0)


    def get_ref_names(self, bam_file):
        """Get the names of reference sequences and reference genomes.

        Parameters
        ----------
            bam_file : str
                Path to bam file.

        Returns
        -------
            ref_seq : set[str]
                A set of reference sequences ids
            ref_genomes : set[str]
                A set of reference genome ids
        """
        ref_seqs = set()
        ref_genomes = set()

        infile = pysam.AlignmentFile(bam_file, "rb")
        for sq in infile.header["SQ"]:
            ref_seq_id = sq["SN"]

            if ref_seq_id not in ref_seqs:
                ref_seqs.add(ref_seq_id)

            # see if we can extract a reference genome id
            match = re.match(self.regex, ref_seq_id)
            if match is None:
                ref_genome = ref_seq_id
            else:
                ref_genome = match.group(0)

            if ref_genome not in ref_genomes:
                ref_genomes.add(ref_genome)

        return ref_seqs, ref_genomes


    def extract_read_pos(self, bam_file):
        """Get the read positions from sam_file for each reference genome.

        Parameters
        ----------
            bam_file : str
                Path to a bam file.

        Returns
        -------
            read_positions : dict[str -> list]
                A dictionary whose key is a reference sequence id,
                and the value is a list of read coordinates along that
                reference sequence
            lengths : dict[str -> int]
                A dictionary whose key is a reference sequence id,
                and the value is the length of that reference sequence
        """

        # reference sequence id to list of read positions
        read_positions = {}
        num_mapped_reads = {}

        # mate pairs both have the same query name,
        # so if we only want to select one read per mate,
        # we need to quickly check duplicate read names
        read_names = set()

        alignment_file = pysam.AlignmentFile(bam_file, "rb")
        for aln in alignment_file:
            if aln.is_unmapped:
                continue

            read_name = aln.query_name
            ref_name = aln.reference_name
            pos = aln.get_reference_positions()[0]
            qual_scores = aln.query_qualities
            mean_qual = np.mean(qual_scores)
            mapping_quality = aln.mapping_quality

            if mean_qual < self.min_avg_qual:
                continue

            # we've already seen this read's mate pair
            if read_name in read_names:
                continue

            if ref_name in read_positions:
                read_positions[ref_name].append(pos)
            else:
                read_positions[ref_name] = [pos]
            num_mapped_reads[ref_name] = num_mapped_reads.get(ref_name, 0) + 1

            read_names.add(read_name)

        lengths = self.get_ref_seq_lengths(bam_file)
        for ref_name in lengths:
            if ref_name not in read_positions:
                read_positions[ref_name] = []
                num_mapped_reads[ref_name] = 0

        return read_positions, lengths


    def get_ref_seq_lengths(self, bam_file):
        """Extract the lengths of each reference sequence from a bam file.

        Parameters
        ----------
            bam_file : str
                Path to a bam file.

        Returns
        -------
            ref_seq_lengths : dict[str -> int]
                A dictionary whose key is a reference sequence id,
                and value is the length of that reference sequence.
        """
        ref_seq_lengths = {}
        alignment_file = pysam.AlignmentFile(bam_file, "rb")
        for sq in alignment_file.header["SQ"]:

            if "LN" not in sq:
                continue
            ref_seq_lengths[sq["SN"]] = sq["LN"]

        return ref_seq_lengths


    def process_bam(self, bam_file):
        """Extract read coordinates along each reference sequence.

        Parameters
        ----------
            bam_file : str
                Path to a bam file.

        Returns
        -------
            coverage_maps : dict[str -> CoverageMap]
                A dictionary whose key is a reference genome id,
                and value is a CoverageMap.
        """
        print_info("BamProcessor", "processing {}".format(bam_file))

        # extract read positions. each variable is a dictionary where the
        # key is the reference sequence name
        read_positions, lengths = self.extract_read_pos(bam_file)

        # we now need to group reference sequencies by species
        # assemblies will have multiple contigs, while complete reference
        # genomes will have a single contig
        contig_read_positions = {}
        contig_lengths = {}

        # group by reference genome
        # for contigs, this will group all contigs together by species
        for ref in sorted(read_positions):
            match = re.match(self.regex, ref)

            if match is not None:
                genome_id = match.group(0)
            else:
                genome_id = ref

            if genome_id not in contig_read_positions:
                contig_read_positions[genome_id] = {}
                contig_lengths[genome_id] = {}

            contig_read_positions[genome_id][ref] = read_positions[ref]
            contig_lengths[genome_id][ref] = lengths[ref]

        # now, store coverage maps for each refence genome
        coverage_maps = {}
        for genome_id in contig_read_positions:
            is_assembly = len(contig_read_positions[genome_id]) > 1

            if is_assembly:
                coverage_maps[genome_id] = \
                    CoverageMapContig(bam_file, genome_id, contig_read_positions[genome_id], contig_lengths[genome_id])
            else:
                ref_id = list(contig_read_positions[genome_id].keys())[0]
                coverage_maps[genome_id] = \
                    CoverageMapRef(bam_file, genome_id, contig_read_positions[genome_id][ref_id], contig_lengths[genome_id][ref_id])

        return coverage_maps


    def merge(self, bam_files, out_bam):
        """Merge many bam files from different indexes into one, taking the
        reads with the highest mapping quality from each bam.

        Parameters
        ----------
            bam_files : list[str]
                A list of bam files to merge.
            out_bam : str
                Location to store the merged bam file.
        """
        print_info("BamProcessor", "merging bam_files {}".format(bam_files))
        print_info("BamProcessor", "finding reads with highest mapq")
        # bam header: SN => LN
        seq_len = {}
        # read_id => (best_bamfile, best_score)
        mapq = {}
        # need to find the bam file with the highest mapping quality
        for bam_file in bam_files:
            bf = pysam.AlignmentFile(bam_file, "rb")

            for sq in bf.header["SQ"]:
                seq_len[sq["SN"]] = sq["LN"]

            for aln in bf:
                if aln.is_unmapped:
                    continue

                if aln.query_name not in mapq:
                    mapq[aln.query_name] = (bam_file, aln.mapping_quality)
                elif mapq[aln.query_name][1] < aln.mapping_quality:
                    mapq[aln.query_name] = (bam_file, aln.mapping_quality)

            bf.close()

        # combined header
        header = { 
            "HD" : {"VN": "1.0", "SO": "unsorted"},
            "SQ" : [{"SN": sq, "LN" : seq_len[sq]} for sq in seq_len]
        }

        print_info("BamProcessor", "writing merged file {}".format(out_bam))
        out = pysam.AlignmentFile(out_bam, "wb", header=header)
        for bam_file in bam_files:
            bf = pysam.AlignmentFile(bam_file, "rb")
            for aln in bf:
                if aln.is_unmapped:
                    continue

                if bam_file == mapq[aln.query_name][0]:
                    out.write(aln)
            bf.close()
        print_info("BamProcessor", "finished writing {}".format(out_bam))




class CoverageMap:
    """Data structure to store read positions along reference genomes.

    Parameters
    ----------
        bam_file : str
            Bam file used to construct coverage map
        genome_id : str
            Unique ID of the reference genome
        is_assembly : bool
            If true, the reference genome is an assembly (and has multiple
            contigs). Otherwise is a complete reference genome with a single
            contig.
    """

    def __init__(self, bam_file, genome_id, is_assembly):
        self.bam_file = bam_file
        self.genome_id = genome_id
        self.sample_id = os.path.basename(bam_file).split(".")[0]
        self.is_assembly = is_assembly



class CoverageMapRef(CoverageMap):
    """Data structure to store read positions from complete reference genomes.

    Parameters
    ----------
        bam_file : str
            Bam file used to construct coverage map
        genome_id : str
            Unique ID of the reference genome
        read_positions : list[int]
            A list of coordinates, one per read.
    """

    def __init__(self, bam_file, genome_id, read_positions, length):
        super().__init__(bam_file, genome_id, is_assembly=False)
        self.read_positions = read_positions
        self.length = length


    def get_length(self):
        """Get the lenth of the genome.
        
        Returns
        -------
            length : int
                The length of the reference genome
        """
        return self.length


    def get_reads(self):
        """Get the read coordinates along the reference genome.

        Returns
        -------
            reads : numpy.array
                A numpy array containing the coordinates of each read
        """
        reads = np.copy(self.read_positions)
        return reads



class CoverageMapContig(CoverageMap):
    """Data structure to store read positions from assemblies.

    Parameters
    ----------
        bam_file : str
            Bam file used to construct coverage map
        genome_id : str
            Unique ID of the reference genome
        contig_read_positions : dict[str -> list]
            A dictionary whose key is a contig id (sequence id),
            and value is a list of coordinates from the reads
            along that contig.
        contig_lengths: dict[str -> int]
            A dictionary whose key is a contig id (sequence id),
            and value is the length of that contig.
    """

    def __init__(self, bam_file, genome_id, contig_read_positions, contig_lengths):
        super().__init__(bam_file, genome_id, is_assembly=True)
        self.contig_ids = np.sort(list(contig_read_positions.keys()))
        self.contig_read_positions = contig_read_positions
        self.contig_lengths = contig_lengths

        self.binned_reads = {}
        # total number of 10Kb bins
        self.total_bins = None
        # fraction with nonzero read count
        self.frac_nonzero = None
        # total reads prefiltering
        self.reads = None
        # flag for qc check
        self.passed_qc_flag = None


    def get_length(self, contig_id):
        """Get the length of a contig.
        
        Parameters
        ----------
            contig_id : str
                Reference sequence id of the contig

        Returns
        -------
            length : int
                The length of the contig
        """
        length = self.contig_lengths[contig_id]
        return length


    def get_reads(self, contig_id):
        """Get the coordinates of all reads for a contig.
        
        Parameters
        ----------
            contig_id : str
                Reference sequence id of the contig
        
        Returns
        -------
            reads : numpy.array
                The coordinates of each read along the contig
        """
        reads = np.copy(self.contig_read_positions[contig_id])
        return reads


    def bin_reads_10Kb(self, contig_id):
        """Count reads in 10Kb windows.

        Parameters
        ----------
            contig_id : str
                Reference sequence id of the contig

        Returns
        -------
            binned_reads : numpy.array
                Total number of reads in each 10Kb window
        """
        bin_size = 10000
        contig_length = self.contig_lengths[contig_id]

        if contig_id not in self.binned_reads:
            binned_counts = []

            nbins = int(math.floor(contig_length / bin_size))
            bin_counts = np.zeros(nbins)

            for r in self.contig_read_positions[contig_id]:
                rbin = int(math.floor(nbins * r / contig_length))
                bin_counts[rbin] += 1

            self.binned_reads[contig_id] = bin_counts

        return self.binned_reads[contig_id]


    def passed_qc(self):
        """Run basic quality control checks. Sets the passed_gc_flag attribute.
        """
        if self.passed_qc_flag is not None:
            return self.passed_qc_flag
        total_bins = 0
        total_reads = 0
        zero_bins = 0
        for cid in self.contig_ids:
            length = self.contig_lengths[cid]
            if length < 11000: continue
            bins = self.bin_reads_10Kb(cid)
            total_bins += bins.size
            total_reads += bins.sum()
            zero_bins += (bins == 0).sum()

        if zero_bins / total_bins > 0.5 and total_bins > 50:
            self.passed_qc_flag = False
        else:
            self.passed_qc_flag = True

        self.frac_nonzero = zero_bins / total_bins
        self.reads  = total_reads
        self.total_bins = total_bins
        return self.passed_qc_flag