from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import os


def write_xyz(selection: str, universe: mda.Universe, pathname: str):
    group = universe.select_atoms('resid ' + selection)

    pathname += '.xyz'
    with mda.Writer(pathname, group.n_atoms) as coords:
        coords.write(group)


# Gets xyz data from saved coords
def read_xyz(xyz_path):
    with open(xyz_path, 'r') as coord_file:
        xyz_data = []
        coord_data = coord_file.readlines()[1:]
        images = 0
        for line in coord_data:
            if line == '\n':
                images += 1

        for line in coord_data:
            if '.' in line:
                xyz_data.append(line)
        return xyz_data

    # Saves coords as a string


def save_sapt_in(coords0: list, coords1: list, memory: int, path: str):
    coord_data = '0 1\n'

    for line in coords0:
        if '|' not in line:
            items = line.split()
            line = ''
            items[0] = items[0][0]
            for item in items:
                item += (' ' + item)
                coord_data += (line + '\n')

        coord_data += '--\n'
        coord_data += '-1 1\n'

    for line in coords1:
        if '|' not in line:
            items = line.split()
            line = ''
            items[0] = items[0][0]
            for item in items:
                line += (' ' + item)
            coord_data += (line + '\n')

    coord_data += ''
    coord_data += 'units angstrom\n ' \
                  '\n' \
                  '}\n' \
                  '\nset\n' \
                  '{basis jun-cc-pVDZ\n' \
                  'scf_type df\n' \
                  'freeze_core\n' \
                  '{\n'

    coord_data += '\n' + 'memory ' + str(memory) + 'GB\n'

    coord_data += "energy('sapt0')"

    with open(path, 'w+') as input_file:
        for line in coord_data:
            input_file.write(line)


def check_inputs(selection: list, start: int, stop: int, step: int, universe: mda.Universe):
    ag_sel = selection[0]
    ag_names = selection[1]
    ag_pair = selection[2]

    # Testing names and selections
    if len(ag_sel) > ag_pair:
        raise InputError('Not all selections are named')
    elif len(ag_sel) < len(ag_names):
        raise InputError('Too many selection names for number of selections')

    for sel in ag_sel:
        try:
            ag = universe.select_atoms(sel)
        except mda.SelectionError:
            raise InputError('Error in selection:', sel)
        except mda.SelectionWarning:
            raise InputError('Warning on selection: ', sel)

    for pair in ag_pair:
        if len(pair) != 2:
            raise InputError('Pairs must be a python list of string with only two items')
        found0 = False
        found1 = False
        for name in ag_names:
            if pair[0] == name:
                found = True
            if found0 is False:
                raise InputError(pair[0], ' in ', pair, 'group_pair_selections is not in defined in atom_group_names')
            if found1 is False:
                raise InputError(pair[1], ' in ', pair, 'group_pair_selections is not in defined in atom_group_names')

        if start >= stop:
            raise InputError('Start is greater than or equal to stop')
        if step >= stop:
            raise InputError('Step is greater than or equal to stop')

        if len(universe.trajectory) < stop:
            raise InputError('Stop exceeds length of trajectory, trajectory is ', len(universe.trajectory), ' frames')

    print('Input Parameters Accepted')


class InputError(Exception):
    pass


class Psi4SAPTGenerator(AnalysisBase):
    def __init__(self, universe: mda.Universe, selections: list, memory: int):
        super(Psi4SAPTGenerator).__init__(universe.trajectory)
        self._unv = universe
        self._sel = selections
        self._mem = memory

    def _prepare(self):
        self.time_list = []
        self.selection_coords = self._sel[0]
        self.selection_names = self._sel[1]
        self.interaction_pairs = self._sel[2]

    def _single_frame(self):
        for ind in range(len(self.selection_coords)):
            write_xyz(self.selection_coords[ind], self._unv, self.selection_names[ind])

        for pair in self.interaction_pairs:
            coords0 = read_xyz(pair[0] + '.xyz')
            coords1 = read_xyz(pair[1] + '.xyz')

            path = 'frame' + self._unv.trajectory.time + '_' + pair[0] + '_' + pair[1]
            save_sapt_in(coords0, coords1, self._mem, path)

    def _conclude(self):
        for path in self.selection_names:
            path += '.xyz'
            os.remove(path)
