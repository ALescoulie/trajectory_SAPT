from inputwriter import *
import MDAnalysis as mda
import unittest


class FrameTest(unittest.TestCase):

    def test_read_write_xyz0(self) -> None:
        top = 'testfiles/testtop.psf'
        trj = 'testfiles/testtraj.dcd'
        unv = mda.Universe(top, trj)
        ag = 'resid 1'
        path = 'Met1'
        write_xyz(ag, unv, path)
        xyz_coords = read_xyz(path + '.xyz')
        with open(path + '.xyz', 'r') as testfile:
            test_data = testfile.readlines()[2:]
            ind = 0
            for line in test_data:
                if '.' in line:
                    self.assertEqual(line, xyz_coords[ind])
                    ind += 1

    def test_read_write_xyz1(self) -> None:
        top = 'testfiles/testtop.psf'
        trj = 'testfiles/testtraj.dcd'
        unv = mda.Universe(top, trj)
        ag = 'resid 2'
        path = 'Arg2'
        write_xyz(ag, unv, path)
        xyz_coords = read_xyz(path + '.xyz')
        with open(path + '.xyz', 'r') as testfile:
            test_data = testfile.readlines()[2:]
            ind = 0
            for line in test_data:
                if '.' in line:
                    self.assertEqual(line, xyz_coords[ind])
                    ind += 1

    def test_save_sapt_in(self) -> None:
        coords1 = read_xyz('Met1.xyz')
        coords2 = read_xyz('Arg2.xyz')
        save_sapt_in(coords1, coords2, 12, 'frame1_Met1_Arg2.in')
        sapt1 = []
        sapt2 = []
        with open('frame1_Met1_Arg2.in', 'r') as sapt:
            sapt_data = sapt.readlines()
            for line in sapt_data:
                sapt1.append(line)

        with open('testfiles/test.in', 'r') as test:
            test_data = test.readlines()
            for line in test_data:
                sapt2.append(line)
        self.assertEqual(sapt1, sapt2)


if __name__ == '__main__':
    unittest.main()
