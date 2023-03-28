from pathlib import Path
from Bio import SeqIO
import numpy as np
from fragment_analyzer.method.method import Method
from fragment_analyzer.utils.baseline_removal import baseline_arPLS


class Fsa:
    """
    A class to represent a Fragment analysis data object.

    ...

    Attributes
    ----------
    file : str
        file name
    method : Method object
        Pydantic method object  to be used
    correct_baseline : bool
        determines whether the baseline should be corrected

    Methods
    -------
    get_channel_names():
        Returns names of channels

    get_channels():
        Returns names of channels

    _set_lookups():
        populates dicts self.name_to_channel and self.channel_to_name

    _apply_method():
        calls self._set_ladder_data and self._set_raw_traces
        optionally calls self._set_baseline_corrected_traces

    _set_raw_traces():
        populates dict self.traces_raw

    set_method():
        sets the dicts self.traces_raw and self.traces_corrected_baseline
        sets the string variables self.ladder_channel self.ladder_name
        sets the string variable self.method_name
        calls the methods self._apply_method and self._set_lookups

    """
    def __init__(
        self,
        file: str,
        method: Method = None,
        correct_baseline: bool = False,
    ) -> None:

        self.file = Path(file)
        self.method = None
        self.correct_baseline = correct_baseline

        if method is not None:
            self.method_name = method.name
            self.method = method

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

    def get_channels_names(self):
        return self.name_to_channel.keys()

    def _set_lookups(self):
        for channel in self.method.channels:
            self.name_to_channel[channel.channel_name] = channel.channel
            self.channel_to_name[channel.channel] = channel.channel_name

    def _apply_method(self):
        self._set_ladder_data()
        self._set_raw_traces()
        if self.correct_baseline is True:
            self._set_baseline_corrected_traces()

    def _set_ladder_data(self):
        for channel in self.method.channels:
            if channel.ladder:
                self.ladder_channel = channel.channel
                self.ladder_name = channel.channel_name

    def _set_baseline_corrected_traces(self):
        for channel in self.method.channels:
            self.traces_corrected_baseline[channel.channel_name] = np.array(
                baseline_arPLS(self.abif_raw[channel.channel])
            )

    def _set_raw_traces(self):
        for channel in self.method.channels:
            self.traces_raw[channel.channel_name] = np.array(
                self.abif_raw[channel.channel]
            )

    def set_method(self, method: dict):
        self.traces_raw = {}
        self.traces_corrected_baseline = {}
        self.ladder_channel = None
        self.ladder_name = None

        self._apply_method()
        self._set_lookups()

    def __repr__(self):
        return f'Fsa(file=\'{self.file}\', filename=\'{self.filename}\', method=\'{self.method}\', method_name=\'{self.method_name}\', ' \
               f'corrected_baseline=\'{self.correct_baseline}\', ladder_channel=\'{self.ladder_channel}\', ladder_name=\'{self.ladder_name}\', ' \
               f'name_to_channel=\'{self.name_to_channel}\', channel_to_name=\'{self.channel_to_name}\')'

