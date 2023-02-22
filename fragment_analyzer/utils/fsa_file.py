from pathlib import Path
from Bio import SeqIO
import numpy as np

from fragment_analyzer.ladders.ladders import LADDERS
from fragment_analyzer.utils.baseline_removal import baseline_arPLS


class FsaFile:
    def __init__(
        self,
        file: str,
        ladder: str,
        normalize: bool = False,
        trace_channel: str = "DATA1",
    ) -> None:
        peak_count_padding = 3
        self.file = Path(file)
        self.file_name = self.file.parts[-1]

        self.fsa = SeqIO.read(file, "abi").annotations["abif_raw"]
        self.ladder = ladder
        self.trace_channel = trace_channel
        self.normalize = normalize

        self.size_standard_channel = LADDERS[ladder]["channel"]
        self.ref_sizes = LADDERS[ladder]["sizes"]
        self.ref_count = self.ref_sizes.size
        self.min_interpeak_distance = LADDERS[ladder]["distance"]
        self.min_height = LADDERS[ladder]["height"]
        self.max_ladder_trace_distance = LADDERS[ladder]["max_ladder_trace_distance"]
        self.max_peak_count = self.ref_sizes.size + peak_count_padding

        if normalize:
            self.size_standard = np.array(
                baseline_arPLS(self.fsa[self.size_standard_channel])
            )
            self.trace = np.array(baseline_arPLS(self.fsa[trace_channel]))
        else:
            self.size_standard = np.array(self.fsa[self.size_standard_channel])
            self.trace = np.array(self.fsa[trace_channel])

    def __repr__(self):
        representer = f"""
            FsaFile-object with following parameters:
            
            File: {self.file}
            Filename: {self.file_name}
            Size standard channel: {self.size_standard_channel}
            Ladder name: {self.ladder}
            Number of ladder steps: {self.ref_count}
            Minimum interpeak distance: {self.min_interpeak_distance}
            Minimum height: {self.min_height}
            Minimum Ladder trace distance: {self.max_ladder_trace_distance}
            Maximum peak count: {self.max_peak_count}
            Normalized data: {self.normalize}
            Trace Channel: {self.trace_channel}
            Ladder Sizes: {self.ref_sizes}
            """
        return representer
