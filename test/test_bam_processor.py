import os
import unittest

import pysam

from src.coptr.bam_processor import BamProcessor


def read_to_dict(read):
    read = read.strip("\n").split()
    read_dict = {
        "name": read[0],
        "flag": read[1],
        "ref_name": read[2],
        "ref_pos": read[3],
        "map_quality": read[4],
        "cigar": read[5],
        "next_ref_name": read[6],
        "next_ref_pos": read[7],
        "length": read[8],
        "seq": read[9],
        "qual": read[10],
        "tags": read[11:],
    }
    return read_dict


class TestBamProcessor(unittest.TestCase):
    def setUp(self):
        header = {
            "HD": {"VN": "1.0", "SO": "unsorted"},
            "SQ": [
                {"SN": "ref1|seq1", "LN": 1000000},
                {"SN": "ref1|seq2", "LN": 1000000},
                {"SN": "ref2|seq1", "LN": 1000000},
            ],
        }

        bam1_reads = [
            "read1     0       ref1|seq1      24975   42      80M     *       0       0       TGGGCCAGAAAAAATGACTTCTCCATCTCGCTGCCGGTAGACCGACTCTCTTTTCTGCTGGCGGTTGCCACGCTGAGCGG        AAAAAF.A.FFAFFFFFAFFFFFFFFFFFFFF<FFFFAFFFFFFA.FFFFA<7FFFFFFFF<FFFFFF))<FFFFF.FFF        AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:76A3       YT:Z:UU",
            "read2     0       ref1|seq1      20984   42      80M     *       0       0       GTTTAAACAGTTGTTGTTGTTCTTCCTGCGATACTCCACTTCCAGAAGCCATAATCGTCATTTTGATAACAGCGTGGTTG        AAAAA.<FFAFFFFFFF<FFAFF)FFFFF<FFF.FFA)FFAF<F<F<.FF<F.FFAFFF7FAFFF.AF.<)F7FFAAFFF        AS:i:-6 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:33A13T32   YT:Z:UU",
            "read3     0       ref2|seq1       3210    42      80M     *       0       0       ACCTACCACTTCACCGACATATTCATGGCCCACGACCATCGGCACCGGGATGGATTTTTGCGACCACTCATCCCAGTTAT        AAAA7FAFFFFF.FFFFF<FFFAA7FFFFFF7FFFFFFFA<FF7FFAF<F.FF.FFF7FFFFAF<FFFFAFFFFA77FFF        AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:53T26      YT:Z:UU",
            "read4     0       ref1|seq2       9298    23      79M     *       0       0       CAGCATCGCTTCCAAAAATAGTAGTGCAGTTGATCGGAGTAGGAGCGTAATGGATTGCCTGCGTGATTGGCTATCTGGC AAAAAF.A.FFAFFFFFAFFFFFFFFFFFFFF<FFFFAFFFFFFA.FFFFA<7FFFFFFFF<FFFFFF))<FFFFF.FF AS:i:-23        XN:i:0  XM:i:6  XO:i:0  XG:i:0  NM:i:6  MD:Z:19T8A0C2T4T10T30   YT:Z:UU",
        ]

        aln_header = pysam.AlignmentHeader().from_dict(header)
        aln_segment = pysam.AlignedSegment()
        test_bam1 = pysam.AlignmentFile("test/test_bam1.bam", "wb", header=header)
        for read in bam1_reads:
            read = read_to_dict(read)
            test_bam1.write(aln_segment.from_dict(read, aln_header))
        test_bam1.close()

        bam2_reads = [
            "read1     0       ref2|seq1      24975   50      80M     *       0       0       TGGGCCAGAAAAAATGACTTCTCCATCTCGCTGCCGGTAGACCGACTCTCTTTTCTGCTGGCGGTTGCCACGCTGAGCGG        AAAAAF.A.FFAFFFFFAFFFFFFFFFFFFFF<FFFFAFFFFFFA.FFFFA<7FFFFFFFF<FFFFFF))<FFFFF.FFF        AS:i:0 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:76A3       YT:Z:UU",
            "read2     0       ref2|seq1      20984   30      80M     *       0       0       GTTTAAACAGTTGTTGTTGTTCTTCCTGCGATACTCCACTTCCAGAAGCCATAATCGTCATTTTGATAACAGCGTGGTTG        AAAAA.<FFAFFFFFFF<FFAFF)FFFFF<FFF.FFA)FFAF<F<F<.FF<F.FFAFFF7FAFFF.AF.<)F7FFAAFFF        AS:i:-12 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:33A13T32   YT:Z:UU",
        ]

        test_bam2 = pysam.AlignmentFile("test/test_bam2.bam", "wb", header=header)
        for read in bam2_reads:
            read = read_to_dict(read)
            test_bam2.write(aln_segment.from_dict(read, aln_header))
        test_bam2.close()

    def tearDown(self):
        os.remove("test/test_bam1.bam")
        os.remove("test/test_bam2.bam")

    def test_process_bam(self):
        bam_processor = BamProcessor()
        coverage_maps = bam_processor.process_bam("test/test_bam1.bam")

        genome1 = "ref1"
        genome2 = "ref2"

        # check that reference genome ids were processed properly
        self.assertTrue(genome1 in coverage_maps)
        self.assertTrue(genome2 in coverage_maps)

        # check that assemblies and compete genomes are identified
        self.assertTrue(coverage_maps[genome1].is_assembly)
        self.assertFalse(coverage_maps[genome2].is_assembly)

        # check that reads are stored
        self.assertTrue(3210 - 1 in coverage_maps[genome2].read_positions)
        self.assertTrue(
            24975 - 1 in coverage_maps[genome1].contig_read_positions["ref1|seq1"]
        )
        self.assertTrue(
            20984 - 1 in coverage_maps[genome1].contig_read_positions["ref1|seq1"]
        )
        self.assertTrue(
            9298 - 1 in coverage_maps[genome1].contig_read_positions["ref1|seq2"]
        )

    def test_merge_bam(self):
        bam_processor = BamProcessor()
        bam_processor.merge(
            ["test/test_bam1.bam", "test/test_bam2.bam"], "test/test_bam_merged.bam"
        )

        infile = pysam.AlignmentFile("test/test_bam_merged.bam", "rb")
        nreads = 0
        for read in infile:
            nreads += 1
            if read.query_name == "read1":
                self.assertTrue(read.reference_name == "ref2|seq1")
            elif read.query_name == "read2":
                self.assertTrue(read.reference_name == "ref1|seq1")
            elif read.query_name == "read3":
                self.assertTrue(read.reference_name == "ref2|seq1")
            elif read.query_name == "read4":
                self.assertTrue(read.reference_name == "ref1|seq2")

        self.assertTrue(nreads == 4)
        os.remove("test/test_bam_merged.bam")


if __name__ == "__main__":
    unittest.main()
