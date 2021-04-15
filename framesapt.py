from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import psi4
import os


def write_xyz(selection, universe):
    group = universe.select_atoms('resid')
    items = selection.split()
    name = ''
    for item in items:
        if item != 'and' or item != 'resid' or item != 'name':
            name += '_' + item
    with mda.Writer((name + '.xyz'), group.n_atoms) as coords:
        coords.write(group)
    return


def read_xyz(selection):
    path = selection + '.xyz'
    xyz_data = []
    with open(path, 'r') as coord_file:
        coord_data = coord_file.readlines()[1:]
        images = 0
        for line in coord_data:
            if line == '\n':
                images += 1

        for line in coord_data:
            if '.' in line:
                xyz_data.append(line)
    return xyz_data


def save_coords(coords0, coords1):
    coord_data = '0 1\n'

    for line in coords0:
        coord_data += line

    coord_data += '--\n'
    coord_data += '0 1\n'

    for line in coords1:
        coord_data += line

    coord_data += 'units angstrom'
    return coord_data


def calculate_spat(selections: list, universe: mda.Universe, frame: int):
    for pair in selections:
        write_xyz(pair[0], universe)
        write_xyz(pair[1], universe)

        group0_coords = read_xyz(pair[0])
        group1_coords = read_xyz(pair[1])

        psi_in = save_coords(group0_coords, group1_coords)
        dimer = psi4.geometry(psi_in)

        psi4.set_options({'scf_type': 'df',
                          'freeze_core': 'true'})
        psi4.energy('spat0/jun-cc-pvdz', molecule=dimer)


class FrameSAPT(AnalysisBase):
    def __init__(self, universe, selections):
        super(FrameSAPT).__init__(universe.trajectory)
        self._unv = universe
        self._sel = selections

    def _prepare(self):
        self.time_list = []

    def _single_frame(self):
        pass

    def _conclude(self):
        pass


if __name__ == '__main__':
    unv = mda.Universe()
    atp = unv.select_atoms('resid 378')
    er2k = unv.select_atoms('segid 2')
    
