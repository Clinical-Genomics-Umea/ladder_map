import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from scipy import stats, signal
from Bio import SeqIO
from sklearn.linear_model import LinearRegression
from pathlib import Path
from sklearn.metrics import r2_score
from sklearn.preprocessing import minmax_scale
from scipy.interpolate import UnivariateSpline

from .ladders.ladders import LADDERS, CHANNELS
from .baseline_removal import baseline_arPLS


class TracePeakLadderAssigner:
    def __init__(self, ladder: str, normalize_peaks: bool = False) -> None:
        
        self.channel = CHANNELS[ladder]
        self.ref_sizes = LADDERS[ladder]['sizes']
        self.ref_count = self.ref_sizes.size
        self.normalize_peaks = normalize_peaks
        self.max_peak_count = self.ref_sizes.size + 3

        self.min_interpeak_distance = LADDERS[ladder]['distance']
        self.height = LADDERS[ladder]['height']
        self.max_ladder_trace_distance = LADDERS[ladder]['max_ladder_trace_distance']

        # self.data = SeqIO.read(data_, "abi").annotations["abif_raw"]

        # if self.normalize_peaks:
        #     self.ladder_trace = np.array(baseline_arPLS(self.data[self.channel]))
        # else:
        #     self.ladder_trace = np.array(self.data[self.channel])

        # self.max_peak_count = max_peak_count
        # self.distance = distance
        # self.height = height
        # self.peaks = self.get_peaks()
        # self.max_diff = np.min(
        #     [np.diff(self.peaks).max() * max_diff_coefficient, 300]
        # )  # max_diff can maximum be 300
        # self.graph = self.generate_graph()
        # self.best_ladder_peak_correlation()
        # self._fit_linear_model()

    def assign_ladder_peak_sizes(self, fsa: Path):
        ladder_trace = SeqIO.read(fsa, "abi").annotations["abif_raw"]

        if self.normalize_peaks:
            ladder_trace = np.array(baseline_arPLS(ladder_trace[self.channel]))
        else:
            ladder_trace = np.array(ladder_trace[self.channel])

        peaks = self.get_peaks(ladder_trace)
        graph = self.generate_graph(peaks)
        combinations = self.generate_combinations(graph)

        score, comb = self.get_best_fit(combinations)
        return score, comb

    def get_best_fit(self, combinations, method='2nd_derivative'):
        df = pd.DataFrame()
        df["combination"] = list(combinations)

        if method == "2nd_derivative":
            df['score'] = np.vectorize(self._max_spline_second_derivative_score)(df['combination'])

            df = df.sort_values(by="score", ascending=True)

            print(df.head().to_string())

            best = df.head(1)

            return best.score.squeeze(), best.combination.squeeze()

    def get_peaks(self, ladder_trace) -> np.array:

        peaks_obj = signal.find_peaks(
            ladder_trace, distance=self.min_interpeak_distance, height=self.height
        )

        heights = peaks_obj[1]["peak_heights"]
        peaks = peaks_obj[0]

        df = pd.DataFrame({"peaks": peaks, "heights": heights})

        idxmax = df["heights"].idxmax()
        df = df.drop(idxmax)

        peaks_adj = df.nlargest(self.max_peak_count, ["heights"])

        return peaks_adj["peaks"].sort_values().to_numpy()

    def generate_graph(self, peaks) -> nx.DiGraph:
        G = nx.DiGraph()

        for p in peaks:
            G.add_node(p)

        i = 0
        while i < peaks.size:
            j = i + 1
            while j < peaks.size:
                diff = peaks[j] - peaks[i]
                if diff <= self.max_ladder_trace_distance:
                    G.add_edge(peaks[i], peaks[j], length=diff)
                j += 1
            i += 1

        return G

    def generate_combinations(self, graph):
        start_nodes = [
            node for node in graph.nodes if graph.in_degree(node) == 0
        ]
        end_nodes = [
            node for node in graph.nodes if graph.out_degree(node) == 0
        ]

        all_paths = list(nx.all_simple_paths(graph, start_nodes[0], end_nodes[0]))

        for p_arr in all_paths:
            for i in range(0, len(p_arr) - self.ref_count + 1):
                yield np.array(p_arr[i: i + self.ref_count])


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

    @staticmethod
    def _polynomial_model_inv_r2_score(ladder: np.array, comb: np.array) -> float:
        fit = np.polyfit(ladder, comb, 3)
        predict = np.poly1d(fit)

        return 1 - r2_score(comb, predict(ladder))

    @staticmethod
    def _max_fractional_deviation_score(ladder: np.ndarray, comb: np.ndarray):
        l_intervals = np.diff(ladder)
        c_intervals = np.diff(minmax_scale(comb, feature_range=(ladder[0], ladder[-1])))

        frac_deviation = np.abs(l_intervals - c_intervals) / l_intervals
        max_frac_deviation_idx = np.argmax(frac_deviation)

        return frac_deviation[max_frac_deviation_idx]

    def _max_first_derivative_score(self, combination: np.ndarray):

        comb_scaled = minmax_scale(combination, feature_range=(self.ref_sizes[0], self.ref_sizes[-1]))

        diff_intervals = np.diff(comb_scaled) - np.diff(self.ref_sizes)
        abs_diff_intervals_gradient = np.abs(np.gradient(diff_intervals))
        max_abs_diff_intervals_gradient_idx = np.argmax(abs_diff_intervals_gradient)

        return abs_diff_intervals_gradient[max_abs_diff_intervals_gradient_idx]

    def _max_second_derivative_score(self, combination: np.ndarray):

        comb_scaled = minmax_scale(combination, feature_range=(self.ref_sizes[0], self.ref_sizes[-1]))

        diff_intervals = np.diff(comb_scaled) - np.diff(self.ref_sizes)
        abs_second_derivative = np.abs(np.gradient(np.gradient(diff_intervals)))
        max_second_derivative_idx = np.argmax(abs_second_derivative)

        return abs_second_derivative[max_second_derivative_idx]

    def _max_spline_second_derivative_score(self, combination: np.ndarray):
        spl = UnivariateSpline(self.ref_sizes, combination, s=0)
        der2 = spl.derivative(n=2)
        return max(abs(der2(self.ref_sizes)))


    @property
    def plot_best_sample_ladder(self):
        fig = plt.figure(figsize=(20, 10))
        plt.plot(self.ladder_trace)
        plt.plot(
            self.best_correlated_peaks,
            self.ladder_trace[self.best_correlated_peaks],
            "o",
        )
        plt.title(f"Correlation with Ladder: {self.best_correlation * 100: .2f}")

        for peak, ladder in zip(self.best_correlated_peaks, self.ladder):
            plt.text(peak, self.ladder_trace[peak], ladder)

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
