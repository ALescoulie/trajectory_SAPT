from trajectory_SAPT.inputwriter import *
import MDAnalysis as mda
import unittest


class FrameTest(unittest.TestCase):

    def test_read_write_xyz0(self) -> None:
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
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
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
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
        save_sapt_in(coords1, coords2, '0 1', '1 1', '12', 'frame1_Met1_Arg2.in', 'test')
        sapt1 = []
        sapt2 = []
        with open('frame1_Met1_Arg2.in', 'r') as sapt:
            sapt_data = sapt.readlines()
            for line in sapt_data:
                sapt1.append(line)

        with open('test.in', 'r') as test:
            test_data = test.readlines()
            for line in test_data:
                sapt2.append(line)
        self.assertEqual(sapt1, sapt2)

    def test_check_inputs0(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1']
        ag_pairs = [['Met1', 'Met1', '0 1', '0 1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 1, 4, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), 'Not all selections are named')

    def test_check_inputs1(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1', 'Met3', 'Arg2']
        ag_pairs = [['Met1', 'Met1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 1, 4, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), 'Too many selection names for number of selections')

    def test_check_inputs2(self) -> None:
        atom_groups = ['resid 1', 'res 2']
        group_names = ['Met1', 'Arg2']
        ag_pairs = [['Met1', 'Arg2']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 1, 4, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), 'Error in selection: ' + str(atom_groups[1]))

    def test_check_inputs3(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1', 'Arg2']
        ag_pairs = [['Met1', 'Arg2', 'Met1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 1, 4, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), 'Pairs must be a python list of string with 4 items')

    def test_check_inputs4(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1', 'Arg2']
        ag_pairs = [['Asp2', 'Met1', '1 1', '0 1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 1, 4, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), f'{ag_pairs[0][0]} in {ag_pairs[0]} group_pair_selections is not in defined in atom_group_names')

    def test_check_inputs5(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1', 'Arg2']
        ag_pairs = [['Arg2', 'Asp1', '1 1', '0 1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 1, 4, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), f'{ag_pairs[0][1]} in {ag_pairs[0]} group_pair_selections is not in defined in atom_group_names')

    def test_check_inputs6(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1', 'Arg2']
        ag_pairs = [['Arg2', 'Met1', '1 1', '0 1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 5, 4, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), 'Start is greater than or equal to stop')

    def test_check_inputs7(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1', 'Arg2']
        ag_pairs = [['Arg2', 'Met1', '1 1', '0 1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 5, 10, 15, unv)
        except InputError as e:
            self.assertEqual(str(e), 'Step is greater than or equal to stop')

    def test_check_inputs8(self) -> None:
        atom_groups = ['resid 1', 'resid 2']
        group_names = ['Met1', 'Arg2']
        ag_pairs = [['Arg2', 'Met1', '1 1', '0 1']]
        top = 'testtop.psf'
        trj = 'testtraj.dcd'
        unv = mda.Universe(top, trj)
        selections = [atom_groups, group_names, ag_pairs]
        try:
            check_inputs(selections, 5, 120, 1, unv)
        except InputError as e:
            self.assertEqual(str(e), f'Stop exceeds length of trajectory, trajectory is {len(unv.trajectory)} frames')


if __name__ == '__main__':
    unittest.main()
