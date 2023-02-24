from pathlib import Path
from Bio import SeqIO
import numpy as np

from fragment_analyzer.utils.baseline_removal import baseline_arPLS


class FsaObj:
    def __init__(
        self,
        file: str,
        method_def: dict,
        normalize: bool = False,
    ) -> None:

        self.file = Path(file)
        self.method_def = method_def
        self.file_name = self.file.parts[-1]
        self.fsa = SeqIO.read(file, "abi").annotations["abif_raw"]
        self.normalized = False
        self.traces = {}
        self.ladder_channel_name = None

        for channel_name in self.method_def:
            if self.method_def[channel_name]['ladder'] is True:
                self.ladder_channel = channel_name

        if normalize:
            for channel_name in self.method_def:
                self.traces[channel_name] = np.array(
                        baseline_arPLS(self.fsa[self.method_def[channel_name]['data_channel']])
                )
            self.normalized = True
        else:
            for channel_name in self.method_def:
                self.traces[channel_name] = np.array(
                        self.fsa[self.method_def[channel_name]['data_channel']]
                )

    def __repr__(self):
        ladder_c = self.method_def[self.ladder_channel_name]
        all_channels = [c for c in self.traces]
        r = f"""
            Fsa-object with following parameters:
            
            File: {self.file}
            Filename: {self.file_name}
            Ladder data channel: {ladder_c}
            All trace channels: {all_channels} 
            """
        return r

