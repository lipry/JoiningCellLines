from tqdm import tqdm

class Maps:
    def __init__(self, n):
        self.labels_map = {}
        self.sequences_map = {}
        self.epigenetic_map = {}
        self.n = n

    def build_maps(self, X_sequences_list, X_epigenetic_list, y_list):
        assert len(X_sequences_list) == len(X_epigenetic_list) == len(y_list)

        for idx in tqdm(range(len(X_sequences_list))):

            assert len(X_sequences_list[idx]) == len(X_epigenetic_list[idx]) == len(y_list[idx])

            for s, e, y in zip(X_sequences_list[idx], X_epigenetic_list[idx], y_list[idx]):
                self._add_label(s[0], y, idx)
                self._add_sequence(s[0], s[1])
                self._add_epigenetic_data(s[0], e, idx)

    def _add_label(self, ch, y, cell_line_idx):
        try:
            self.labels_map[ch][cell_line_idx] = y
        except KeyError:
            self.labels_map[ch] = [y if i == cell_line_idx else None for i in range(self.n)]

    def _add_sequence(self, ch, sequence):
        self.sequences_map[ch] = sequence

    def _add_epigenetic_data(self, ch, data, cell_line_idx):
        try:
            self.epigenetic_map[ch][cell_line_idx] = data
        except KeyError:
            self.epigenetic_map[ch] = [data if i == cell_line_idx else None for i in range(self.n)]

    def get_label(self, ch):
        try:
            return self.labels_map[ch]
        except KeyError:
            print("{} not found".format(ch))

    def get_sequence(self, ch):
        try:
            return self.sequences_map[ch]
        except KeyError:
            print("{} not found".format(ch))

    def get_epigenetic_data(self, ch):
        try:
            return self.epigenetic_map[ch]
        except KeyError:
            print("{} not found".format(ch))

