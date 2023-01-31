import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from scipy import stats, signal
from Bio import SeqIO
from sklearn.linear_model import LinearRegression
from scipy.signal import find_peaks, peak_widths
from scipy.integrate import simpson
from lmfit.models import VoigtModel, GaussianModel
import lmfit
from collections import namedtuple
from typing import NamedTuple

from .ladders.ladders import LIZ


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
        self.ladder = LIZ
        self.sample_ladder = np.array(self.data["DATA205"])
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

    def plot_best_sample_ladder(self):
        fig = plt.figure(figsize=(15, 10))
        plt.plot(self.sample_ladder)
        plt.plot(
            self.best_correlated_peaks,
            self.sample_ladder[self.best_correlated_peaks],
            "o",
        )
        plt.title(f"Correlation with Ladder: {self.best_correlation * 100: .2f}")

        for peak, ladder in zip(self.best_correlated_peaks, self.ladder):
            plt.text(peak, self.sample_ladder[peak], ladder)
        return fig


class PeakArea:
    def __init__(
        self, data: pd.DataFrame, start: int, end: int, rel_height: float = 0.95
    ) -> None:
        self.data = data
        self.find_fragment_peaks(start=start, end=end)
        self._find_fragment_peak_width(rel_height=rel_height)
        self.divide_peaks()

    def find_fragment_peaks(
        self, start, end, num_peaks: int = 2, threshold: int = 10000
    ):
        data = self.data.loc[lambda x: x.step_adjusted >= start].loc[
            lambda x: x.step_adjusted <= end
        ]

        peaks, _ = find_peaks(data.peaks, height=threshold)

        # break after 1000 tries
        counter = 0
        while len(peaks) != num_peaks:
            # more peaks
            if len(peaks) > num_peaks:
                threshold += threshold / 2
                peaks, _ = find_peaks(data.peaks, height=threshold)
            # less peaks
            else:
                threshold -= threshold / 2
                peaks, _ = find_peaks(data.peaks, height=threshold)

            counter += 1
            if counter > 1000:
                raise ValueError(
                    "No peaks found in this area. Try another part of the graph."
                )

        self.peaks_index = peaks
        self.peaks_dataframe = data

    def _find_fragment_peak_width(self, rel_height: float = 0.95):
        widths = peak_widths(
            self.peaks_dataframe.peaks,
            self.peaks_index,
            rel_height=rel_height,
        )

        self.left_peak_start, self.right_peak_start = np.floor(widths[2]).astype(int)
        self.left_peak_end, self.right_peak_end = np.ceil(widths[3]).astype(int)
        self.left_peak_height, self.right_peak_height = widths[1]

    def plot_peak_widths(self):
        fig = plt.figure()
        plt.plot(self.peaks_dataframe.step_adjusted, self.peaks_dataframe.peaks)
        plt.plot(
            self.peaks_dataframe.step_adjusted.iloc[self.peaks_index],
            self.peaks_dataframe.peaks.iloc[self.peaks_index],
            "o",
        )

        plt.hlines(
            y=self.left_peak_height,
            xmin=self.peaks_dataframe.step_adjusted.iloc[self.left_peak_start],
            xmax=self.peaks_dataframe.step_adjusted.iloc[self.left_peak_end],
            color="C3",
        )

        plt.hlines(
            y=self.right_peak_height,
            xmin=self.peaks_dataframe.step_adjusted.iloc[self.right_peak_start],
            xmax=self.peaks_dataframe.step_adjusted.iloc[self.right_peak_end],
            color="C4",
        )
        return fig

    def divide_peaks(self):
        # add some padding (10) to the left and right to be sure to include everything in the peak
        self.left_peak_dataframe = self.peaks_dataframe.iloc[
            self.left_peak_start - 10 : self.left_peak_end + 10
        ]
        self.right_peak_dataframe = self.peaks_dataframe.iloc[
            self.right_peak_start - 10 : self.right_peak_end + 10
        ]

    def fit_lmfit_model(self, dataframe: pd.DataFrame, model: lmfit.models):
        dataframe = dataframe.copy()
        y = dataframe.peaks.to_numpy()
        x = dataframe.step_adjusted.to_numpy()

        model = model()
        params = model.guess(y, x)
        out = model.fit(y, params, x=x)

        return dataframe.assign(fitted=out.best_fit), out.values, out.fit_report()

    def plot_lmfit_model(self, model: str):

        if model == "gauss":
            model = GaussianModel
        elif model == "voigt":
            model = VoigtModel
        else:
            raise NotImplementedError(
                "The model is not implemented! Choose between voigt and gauss"
            )

        left_df, left_info, _ = self.fit_lmfit_model(self.left_peak_dataframe, model)
        right_df, right_info, _ = self.fit_lmfit_model(self.right_peak_dataframe, model)

        fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True)
        ax1.plot(left_df.step_adjusted, left_df.peaks, "o")
        ax1.plot(left_df.step_adjusted, left_df.fitted)
        ax1.set_title(f"Left Area: {left_info['amplitude']: .1f}")

        ax2.plot(right_df.step_adjusted, right_df.peaks, "o")
        ax2.plot(right_df.step_adjusted, right_df.fitted)
        ax2.set_title(f"Right Area: {right_info['amplitude']: .1f}")

        fig.suptitle(
            f"Quotient: {right_info['amplitude'] / left_info['amplitude']: .2f}"
        )
        fig.legend(["Raw data", "Model"], bbox_to_anchor=(1.1, 0.9))
        return fig

    def calculate_peak_area(self, function: str = "simpson") -> NamedTuple:
        data = self.peaks_dataframe.peaks.to_numpy()

        if function == "simpson":
            left_area = simpson(data[self.left_peak_start : self.left_peak_end])
            right_area = simpson(data[self.right_peak_start : self.right_peak_end])
            quotient = right_area / left_area
        elif function == "trapz":
            left_area = np.trapz(data[self.left_peak_start : self.left_peak_end])
            right_area = np.trapz(data[self.right_peak_start : self.right_peak_end])
            quotient = right_area / left_area

        Area = namedtuple("Area", "left_area right_area area_quotient")

        return Area(left_area, right_area, quotient)
