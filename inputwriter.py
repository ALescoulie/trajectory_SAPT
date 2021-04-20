def write_xyz(selection, universe):
    group = universe.select_atoms('resid ' + selection)
    path = 'resid' + selection + '.xyz'

    with mda.Writer(path, group.n_atoms) as coords:
        coords.write(group)
    return path


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


def save_sapt_in(coords0, coords1, memory, path):
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
            item = ''
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
    return path


class Psi4SAPTGenerator(AnalysisBase):
    def __init__(self, universe: mda.Universe, selections: list, memory: int):
        super(Psi4SAPTGenerator).__init__(universe.trajectory)
        self._unv = universe
        self._sel = selections

    def _prepare(self):
        self.time_list = []
        self.selection_coords = self._sel[1]
        self.interaction_pairs = self._sel[2]
        self.xyz_paths = []

    def _single_frame(self):
        for group in self.selection_coords:
            self.xyz_paths.append(write_xyz(group, self._unv))

    def _conclude(self):
        pass
