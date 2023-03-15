from pathlib import Path
from Bio import SeqIO
import numpy as np

from fragment_analyzer.utils.baseline_removal import baseline_arPLS


class Fsa:
    def __init__(
        self,
        file: str,
        method: dict = None,
        correct_baseline: bool = False,
    ) -> None:

        self.file = Path(file)
        self.method = None
        self.correct_baseline = correct_baseline

        if method is not None:
            self.method_name = method.keys()[0]
            self.method = method[self.method_name]

        self.filename = self.file.parts[-1]
        self.abif_raw = SeqIO.read(file, 'abi').annotations['abif_raw']

        self.traces_raw = {}
        self.traces_corrected_baseline = {}

        self.ladder_channel = None
        self.ladder_name = None
        self.name_to_channel = {}
        self.channel_to_name = {}

        if method is not None:
            self.method = method
            self._apply_method()
            self._set_lookups()

    def get_channels(self):
        return self.channel_to_name.keys()

    def get_names(self):
        return self.name_to_channel.keys()

    def _set_lookups(self):
        for item in self.method:
            self.name_to_channel[item['name']] = item['channel']
            self.channel_to_name[item['channel']] = item['name']

    def _apply_method(self):
        self._set_ladder_data()
        self._set_raw_traces()
        if self.correct_baseline is True:
            self._set_baseline_corrected_traces()

    def _set_ladder_data(self):
        for item in self.method:
            if item['ladder']:
                self.ladder_channel = item['channel']
                self.ladder_name = item['name']

    def _set_baseline_corrected_traces(self):
        for item in self.method:
            self.traces_corrected_baseline[item['name']] = np.array(
                baseline_arPLS(self.abif_raw[item['channel']])
            )

    def _set_raw_traces(self):
        for item in self.method:
            self.traces_raw[item['name']] = np.array(
                self.abif_raw[item['channel']]
            )

    def set_method(self, method: dict):
        self.traces_raw = {}
        self.traces_corrected_baseline = {}
        self.ladder_channel = None
        self.ladder_name = None

        self.method_name = method.keys()[0]
        self.method = method[self.method_name]
        self._apply_method()
        self._set_lookups()

    def __repr__(self):
        return f'Fsa(file=\'{self.file}\', filename=\'{self.filename}\', method=\'{self.method}\', method_name=\'{self.method_name}\', ' \
               f'corrected_baseline=\'{self.correct_baseline}\', ladder_channel=\'{self.ladder_channel}\', ladder_name=\'{self.ladder_name}\', ' \
               f'name_to_channel=\'{self.name_to_channel}\', channel_to_name=\'{self.channel_to_name}\')'



    #
    # def __repr__(self):
    #     ladder_c = self.method_data[self.ladder_channel]
    #     all_channels = [c for c in self.traces]
    #     r = f"""
    #         Fsa-object with following parameters:
    #
    #         File: {self.file}
    #         Filename: {self.file_name}
    #         Ladder data channel: {ladder_c}
    #         All trace channels: {all_channels}
    #         """
    #     return r

