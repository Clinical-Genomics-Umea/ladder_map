import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from scipy import stats, signal
import pandas as pd
from Bio import SeqIO
from sklearn.linear_model import LinearRegression
from scipy.signal import find_peaks, peak_widths
from scipy.integrate import simpson

from collections import namedtuple
from typing import NamedTuple


class LadderMap:
    def __init__(
        self,
        data: str,
        max_peak_count: int = 38,
        distance: int = 30,
        height: int = 100,
        max_diff_coefficient: float = 1.5,
    ) -> None:
        self.data = SeqIO.read(data, "abi").annotations["abif_raw"]
        self.ladder = np.array(
            [
                100,
                114,
                120,
                140,
                160,
                180,
                200,
                214,
                220,
                240,
                260,
                280,
                300,
                314,
                320,
                340,
                360,
                380,
                400,
                414,
                420,
                440,
                460,
                480,
                500,
                514,
                520,
                540,
                560,
                580,
            ]
        )
        self.sample_ladder = np.array(self.data["DATA205"])
        self.max_peak_count = max_peak_count
        self.distance = distance
        self.height = height
        self.peaks = self.get_peaks()
        self.max_diff = np.diff(self.peaks).max() * max_diff_coefficient
        self.graph = self.generate_graph()
        self.best_ladder_peak_correlation()
        self._fit_linear_model()

    def get_peaks(self) -> np.array:
        peaks_obj = signal.find_peaks(
            self.sample_ladder, distance=self.distance, height=self.height
        )

        heights = peaks_obj[1]["peak_heights"]
        peaks = peaks_obj[0]

        df = pd.DataFrame({"peaks": peaks, "heights": heights})
        idxmax = df["heights"].idxmax()

        idxs_remove = list(range(idxmax + 1))

        df = df.drop(idxs_remove)

        peaks_adj = df.nlargest(self.max_peak_count, ["heights"])

        return peaks_adj["peaks"].sort_values().to_numpy()

    def generate_graph(self) -> nx.DiGraph:
        G = nx.DiGraph()

        for p in self.peaks:
            G.add_node(p)

        i = 0
        while i < self.peaks.size:
            j = i + 1
            while j < self.peaks.size:
                diff = self.peaks[j] - self.peaks[i]
                if diff <= self.max_diff:
                    G.add_edge(self.peaks[i], self.peaks[j], length=diff)
                j += 1
            i += 1

        return G

    def generate_combinations(self):
        start_nodes = [
            node for node in self.graph.nodes if self.graph.in_degree(node) == 0
        ]
        end_nodes = [
            node for node in self.graph.nodes if self.graph.out_degree(node) == 0
        ]

        if len(start_nodes) > 1:
            raise Exception("Can't generate. Too many start nodes.")

        if len(end_nodes) > 1:
            raise Exception("Can't generate, Too many end nodes.")

        all_paths = list(nx.all_simple_paths(self.graph, start_nodes[0], end_nodes[0]))

        for p_arr in all_paths:
            for i in range(0, len(p_arr) - self.ladder.size + 1):
                yield np.array(p_arr[i : i + self.ladder.size])

    def best_ladder_peak_correlation(self):
        result = []
        for combination in self.generate_combinations():
            corr_peaks = stats.pearsonr(self.ladder, combination)

            obj = {"corr_peaks": corr_peaks.statistic, "peaks": combination}

            result.append(obj)

        df = pd.DataFrame(result)
        df = df.sort_values(by="corr_peaks", ascending=False)
        best = df.head(1)

        self.best_correlated_peaks = best.peaks.squeeze()
        self.best_correlation = best.corr_peaks.squeeze()

    def _fit_linear_model(self):
        self.linear_model = LinearRegression()
        self.linear_model.fit(self.best_correlated_peaks.reshape(-1, 1), self.ladder)

    def adjusted_step_dataframe(self, channel: str = "DATA1") -> pd.DataFrame:
        df = (
            pd.DataFrame({"peaks": self.data[channel]})
            .reset_index()
            .rename(columns={"index": "step_raw"})
            .assign(
                step_adjusted=lambda x: self.linear_model.predict(
                    x.step_raw.to_numpy().reshape(-1, 1)
                )
            )
            .loc[lambda x: x.step_adjusted >= 0]
        )

        return df


class PeakArea:
    def __init__(
        self, data: pd.DataFrame, start: int, end: int, rel_height: float = 0.95
    ) -> None:
        self.data = data
        self.find_fragment_peaks(start=start, end=end)
        self.widths = self._find_fragment_peak_width(rel_height=rel_height)

    def find_fragment_peaks(
        self, start, end, num_peaks: int = 2, threshold: int = 10000
    ):
        data = self.data.loc[lambda x: x.step_adjusted >= start].loc[
            lambda x: x.step_adjusted <= end
        ]

        peaks, _ = find_peaks(data.peaks, height=threshold)
        while len(peaks) != num_peaks:
            # more peaks
            if len(peaks) > num_peaks:
                threshold += threshold / 2
                peaks, _ = find_peaks(data.peaks, height=threshold)
            # less peaks
            else:
                threshold -= threshold / 2
                peaks, _ = find_peaks(data.peaks, height=threshold)

        self.peaks_index = peaks
        self.peaks_dataframe = data

    def _find_fragment_peak_width(self, rel_height: float = 0.95):
        widths = peak_widths(
            self.peaks_dataframe.peaks,
            self.peaks_index,
            rel_height=rel_height,
        )
        return widths

    def plot_peak_widths(self):
        self.peaks_dataframe.plot("step_adjusted", "peaks")
        plt.plot(
            self.peaks_dataframe.step_adjusted.iloc[self.peaks_index],
            self.peaks_dataframe.peaks.iloc[self.peaks_index],
            "o",
        )

        plt.hlines(
            y=self.widths[1][0],
            xmin=self.peaks_dataframe.step_adjusted.iloc[
                self.widths[2][0].round().astype(int)
            ],
            xmax=self.peaks_dataframe.step_adjusted.iloc[
                self.widths[3][0].round().astype(int)
            ],
            color="C3",
        )

        plt.hlines(
            y=self.widths[1][1],
            xmin=self.peaks_dataframe.step_adjusted.iloc[
                self.widths[2][1].round().astype(int)
            ],
            xmax=self.peaks_dataframe.step_adjusted.iloc[
                self.widths[3][1].round().astype(int)
            ],
            color="C4",
        )

    def calculate_peak_area(self, function: callable = "simpson") -> NamedTuple:
        data = self.peaks_dataframe.peaks.to_numpy()
        left_start, right_start = self.widths[2].round().astype(int)
        left_end, right_end = self.widths[3].round().astype(int)

        if function == "simpson":
            left_area = simpson(data[left_start:left_end])
            right_area = simpson(data[right_start:right_end])
            quotient = right_area / left_area
        elif function == "trapz":
            left_area = np.trapz(data[left_start:left_end])
            right_area = np.trapz(data[right_start:right_end])
            quotient = right_area / left_area

        Area = namedtuple("Area", "left_area right_area area_quotient")

        return Area(left_area, right_area, quotient)
