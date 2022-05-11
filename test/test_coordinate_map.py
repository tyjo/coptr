import numpy as np
import unittest

from coptr.coptr_ref import CoordinateMap
from coptr.coptr_ref import ReadFilterRef


class TestCoordinateMap(unittest.TestCase):

    def test_coordinate_map1(self):
        rf = ReadFilterRef(0, 0)

        coord = 3
        read_positions = np.array([coord])
        genome_length = 10
        regions = [(0, 1)]

        coord_map = CoordinateMap()
        read_positions, new_genome_length = rf.remove_reads_by_region(read_positions, genome_length, regions, coord_map)

        coord_map = CoordinateMap()
        coord_map.update_break_points(regions[0][1], regions[0][1] - regions[0][0])

        unfiltered_coord = coord_map.translate(read_positions[0])

        self.assertTrue(unfiltered_coord == coord)

    def test_coordinate_map2(self):

        rf = ReadFilterRef(0, 1)

        coordinates = [i for i in range(6, 25)] + \
                      [i for i in range(28, 30)] + \
                      [i for i in range(46, 60)] + \
                      [i for i in range(68, 83)] + \
                      [i for i in range(90, 100)]
        read_positions = np.array(coordinates)
        genome_length = 100
        regions = [(1, 5), (25, 27), (30, 45), (60, 67), (83, 89)]

        coord_map = CoordinateMap()
        new_read_positions, new_genome_length = rf.remove_reads_by_region(read_positions, genome_length, regions, coord_map)
        self.assertTrue(new_read_positions.size == read_positions.size)

        for i,filtered_coord in enumerate(new_read_positions):
            unfiltered_coord = coord_map.translate(filtered_coord)
            self.assertTrue(unfiltered_coord == coordinates[i])



if __name__ == "__main__":
    unittest.main()
