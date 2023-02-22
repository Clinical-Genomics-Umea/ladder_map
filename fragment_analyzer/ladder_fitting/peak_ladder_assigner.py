import pandas as pd
import numpy as np
import networkx as nx
from scipy import signal
from sklearn.metrics import r2_score
from sklearn.preprocessing import minmax_scale
from scipy.interpolate import UnivariateSpline

from fragment_analyzer.utils.fsa_file import FsaFile


class PeakLadderAssigner:
    def __init__(self, fsa_file=FsaFile) -> None:

        self.fsa_file = fsa_file

    def assign_ladder_peak_sizes(self):
        peaks = self.get_peaks(self.fsa_file.size_standard)
        graph = self.generate_graph(peaks)
        combinations = self.generate_combinations(graph)
        best_combination = self.get_best_fit(combinations)

        return best_combination

    def get_peaks(self, size_standard) -> np.array:

        peaks_obj = signal.find_peaks(
            size_standard,
            distance=self.fsa_file.min_interpeak_distance,
            height=self.fsa_file.min_height,
        )

        peaks = peaks_obj[0]
        heights = peaks_obj[1]["peak_heights"]

        df = pd.DataFrame({"peaks": peaks, "heights": heights})

        idxmax = df["heights"].idxmax()
        df = df.drop(idxmax)

        peaks_adj = df.nlargest(self.fsa_file.max_peak_count, ["heights"])

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
                if diff <= self.fsa_file.max_ladder_trace_distance:
                    G.add_edge(peaks[i], peaks[j], length=diff)
                j += 1
            i += 1

        return G

    def generate_combinations(self, graph):
        start_nodes = [node for node in graph.nodes if graph.in_degree(node) == 0]
        end_nodes = [node for node in graph.nodes if graph.out_degree(node) == 0]

        all_paths = list(nx.all_simple_paths(graph, start_nodes[0], end_nodes[0]))

        for p_arr in all_paths:
            for i in range(0, len(p_arr) - self.fsa_file.ref_count + 1):
                yield np.array(p_arr[i : i + self.fsa_file.ref_count])

    def get_best_fit(self, combinations, method="2nd_derivative"):
        df = pd.DataFrame()

        df["combination"] = list(combinations)

        if method == "2nd_derivative":
            df["score"] = np.vectorize(self._max_spline_second_derivative_score)(
                df["combination"]
            )

            df = df.sort_values(by="score", ascending=True)

            best = df.head(1)

            return best.combination.squeeze()

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

        comb_scaled = minmax_scale(
            combination,
            feature_range=(self.fsa_file.ref_sizes[0], self.fsa_file.ref_sizes[-1]),
        )

        diff_intervals = np.diff(comb_scaled) - np.diff(self.fsa_file.ref_sizes)
        abs_diff_intervals_gradient = np.abs(np.gradient(diff_intervals))
        max_abs_diff_intervals_gradient_idx = np.argmax(abs_diff_intervals_gradient)

        return abs_diff_intervals_gradient[max_abs_diff_intervals_gradient_idx]

    def _max_second_derivative_score(self, combination: np.ndarray):

        comb_scaled = minmax_scale(
            combination,
            feature_range=(self.fsa_file.ref_sizes[0], self.fsa_file.ref_sizes[-1]),
        )

        diff_intervals = np.diff(comb_scaled) - np.diff(self.fsa_file.ref_sizes)
        abs_second_derivative = np.abs(np.gradient(np.gradient(diff_intervals)))
        max_second_derivative_idx = np.argmax(abs_second_derivative)

        return abs_second_derivative[max_second_derivative_idx]

    def _max_spline_second_derivative_score(self, combination: np.ndarray):
        spl = UnivariateSpline(self.fsa_file.ref_sizes, combination, s=0)
        der2 = spl.derivative(n=2)
        return max(abs(der2(self.fsa_file.ref_sizes)))
