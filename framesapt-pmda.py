from pmda.parallel import ParallelAnalysisBase
import MDAnalysis as mda
import psi4


def calculate_sapt(selections: list, universe: mda.Universe, frame_num: int, memory: str):
    def write_xyz(selection: str, unv: mda.Universe, frame: int):
        group = universe.select_atoms('resid ' + selection)
        path =  selection + '_' + str(frame) + '.xyz'

        with mda.Writer(path, group.n_atoms) as coords:
            coords.write(group)
        return path

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

    for pair in selections:
        write_xyz(pair[0], universe)
        write_xyz(pair[1], universe)

        group0_coords = read_xyz(pair[0])
        group1_coords = read_xyz(pair[1])

        psi_in = save_coords(group0_coords, group1_coords)
        dimer = psi4.geometry(psi_in)

        psi4.set_options({'scf_type': 'df',
                          'freeze_core': 'true'})
        psi4.energy('sapt0/jun-cc-pvdz', molecule=dimer)


class FrameSAPT(ParallelAnalysisBase):
    def __init__(self, universe, selections):
        super(FrameSAPT).__init__(universe.trajectory)
        self._unv = universe
        self._sel = selections

    def _prepare(self):
        self.time_list = []
        self.energies = []

    def _single_frame(self):
        self.time_list.append()

    def _conclude(self):
        pass


if __name__ == '__main__':
    unv = mda.Universe()
    atp = unv.select_atoms('resid 378')
    er2k = unv.select_atoms('segid 2')
    
