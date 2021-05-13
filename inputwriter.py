from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import os
import sys


def write_xyz(selection: str, universe: mda.Universe, pathname: str):
    group = universe.select_atoms(selection)

    pathname += '.xyz'
    with mda.Writer(pathname, group.n_atoms) as coords:
        coords.write(group)


# Gets xyz data from saved coords
def read_xyz(xyz_path):
    with open(xyz_path, 'r') as coord_file:
        xyz_data = []
        coord_data = coord_file.readlines()[2:]

        for line in coord_data:
            if '.' in line:
                xyz_data.append(line)
        return xyz_data

    # Saves coords as a string


def save_sapt_in(coords0: list, coords1: list, memory: int, path: str, molecule_name: str):
    coord_data = 'molecule %s {\n' % molecule_name

    for line0 in coords0:
        items = line0.split()
        line0 = items[0][0]
        for item in items[1:]:
            line0 += ('\t' + item)
        coord_data += (line0 + '\n')

    coord_data += '--\n'

    for line1 in coords1:
        items = line1.split()
        line1 = items[0][0]
        for item in items[1:]:
            line1 += ('\t' + item)
        coord_data += (line1 + '\n')

    coord_data += '\nunits angstrom\n' \
                  '\n' \
                  '}\n' \
                  '\nset {\n' \
                  'basis jun-cc-pVDZ\n' \
                  'scf_type df\n' \
                  'freeze_core true\n' \
                  '}\n'

    coord_data += '\n' + 'memory ' + str(memory) + 'GB\n'

    coord_data += "\nenergy('sapt0')\n"

    with open(path, 'w+') as input_file:
        for line in coord_data:
            input_file.write(line)


def check_inputs(selection: list, start: int, stop: int, step: int, universe: mda.Universe):
    ag_sel = selection[0]
    ag_names = selection[1]
    ag_pair = selection[2]

    # Testing names and selections
    if len(ag_sel) > len(ag_names):
        raise InputError('Not all selections are named')
    elif len(ag_sel) < len(ag_names):
        raise InputError('Too many selection names for number of selections')

    for sel in ag_sel:
        try:
            ag = universe.select_atoms(sel)
        except mda.SelectionError:
            raise InputError('Error in selection: {}'.format(sel))

    for pair in ag_pair:
        if len(pair) != 2:
            raise InputError('Pairs must be a python list of string with only two items')
        found0 = False
        found1 = False
        for name in ag_names:
            if pair[0] == name:
                found0 = True
            if pair[1] == name:
                found1 = True
        if found0 is False:
            raise InputError(f'{pair[0]} in {pair} group_pair_selections is not in defined in atom_group_names')
        if found1 is False:
            raise InputError(f'{pair[1]} in {pair} group_pair_selections is not in defined in atom_group_names')

        if start >= stop:
            raise InputError('Start is greater than or equal to stop')
        if step >= stop:
            raise InputError('Step is greater than or equal to stop')
        if step == 0:
            raise InputError('Step cannot be 0')

        if len(universe.trajectory) < stop:
            raise InputError(f'Stop exceeds length of trajectory, trajectory is {len(universe.trajectory)} frames')

    print('Input Parameters Accepted')


class InputError(Exception):
    pass


class Psi4SAPTGenerator(AnalysisBase):
    def __init__(self, universe: mda.Universe, selections: list, memory: int, input_directory: str, molecule_name: str):
        super(Psi4SAPTGenerator, self).__init__(universe.trajectory)
        self._unv = universe
        self._sel = selections
        self._mem = memory
        self._dir = input_directory
        self._mol = molecule_name

    def _prepare(self):
        """Defining data structures, selection_coords contains MDAnalysis selection commands,
         selections names contains names of the atom groups selected, and interaction_pairs
         contains the atom group names in list pairs. selections are pre-verified by check_inputs"""
        self.selection_coords = self._sel[0]
        self.selection_names = self._sel[1]
        self.interaction_pairs = self._sel[2]

    def _single_frame(self):
        for ind in range(len(self.selection_coords)):
            write_xyz(self.selection_coords[ind], self._unv, f'{self._dir}/{self.selection_names[ind]}')

        time = int(self._unv.trajectory.time)
        name = f'{self._mol}_{time}'
        for pair in self.interaction_pairs:
            coords0 = read_xyz(f'{self._dir}/{pair[0]}.xyz')
            coords1 = read_xyz(f'{self._dir}/{pair[1]}.xyz')

            path = f'{self._dir}/frame{time}_{pair[0]}_{pair[1]}.in'
            save_sapt_in(coords0, coords1, self._mem, path, name)

    def _conclude(self):
        for path in self.selection_names:
            path = f'{self._dir}/{path}.xyz'
            os.remove(path)


if __name__ == '__main__':
    topology_path = ''
    trajectory_path = ''

