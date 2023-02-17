import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from scipy import stats, signal
from Bio import SeqIO
from sklearn.linear_model import LinearRegression
from pathlib import Path

from .ladders.ladders import LADDERS, CHANNELS
from .baseline_removal import baseline_arPLS


class LadderMap:
    def __init__(
        self,
        data_: str,
        ladder: str,
        normalize_peaks: bool = False,
        max_peak_count: int = 38,
        distance: int = 30,
        height: int = 100,
        max_diff_coefficient: float = 1.5,
    ) -> None:
        self.channel = CHANNELS[ladder]
        self.data_ = Path(data_)
        self.data = SeqIO.read(data_, "abi").annotations["abif_raw"]
        self.ladder = LADDERS[ladder]
        self.normalize_peaks = normalize_peaks

        if self.normalize_peaks:
            self.sample_ladder = np.array(baseline_arPLS(self.data[self.channel]))
        else:
            self.sample_ladder = np.array(self.data[self.channel])

        self.max_peak_count = max_peak_count
        self.distance = distance
        self.height = height
        self.peaks = self.get_peaks()
        self.max_diff = np.min(
            [np.diff(self.peaks).max() * max_diff_coefficient, 300]
        )  # max_diff can maximum be 300
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

        # comment out the below code to make the algorithm work for every file... TODO later
        # if len(start_nodes) > 1:
        #    raise Exception("Can't generate. Too many start nodes.")
        #
        # if len(end_nodes) > 1:
        #    raise Exception("Can't generate, Too many end nodes.")

        # debug
        # seems like this steps takes a long time for some samples...
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

        self.correlation_dataframe = best.explode("peaks").assign(ladder=self.ladder)

    def _fit_linear_model(self):
        self.linear_model = LinearRegression()
        self.linear_model.fit(self.best_correlated_peaks.reshape(-1, 1), self.ladder)

    def adjusted_step_dataframe(self, channel: str = "DATA1") -> pd.DataFrame:
        if self.normalize_peaks:
            data = baseline_arPLS(self.data[channel])
        else:
            data = self.data[channel]

        df = (
            pd.DataFrame({"peaks": data})
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

    @property
    def plot_best_sample_ladder(self):
        fig = plt.figure(figsize=(20, 10))
        plt.plot(self.sample_ladder)
        plt.plot(
            self.best_correlated_peaks,
            self.sample_ladder[self.best_correlated_peaks],
            "o",
        )
        plt.title(f"Correlation with Ladder: {self.best_correlation * 100: .2f}")

        for peak, ladder in zip(self.best_correlated_peaks, self.ladder):
            plt.text(peak, self.sample_ladder[peak], ladder)

        plt.ylabel("intensity")
        plt.xlabel("time")
        plt.grid()
        return fig

    @property
    def plot_ladder_correlation(self):
        fig = plt.figure(figsize=(20, 10))
        plt.plot(
            self.correlation_dataframe.ladder, self.correlation_dataframe.peaks, "o"
        )
        plt.xlabel("Basepairs")
        plt.ylabel("Time")
        plt.grid()
        plt.title("Correlation of found peaks with size-standard")

        return fig
