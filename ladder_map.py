import numpy as np
import networkx as nx
from scipy import stats, signal
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from pandas import DataFrame


def get_peaks(data, max_peak_count=38):
    peaks_obj = signal.find_peaks(data, distance=30, height=100)

    heights = peaks_obj[1]['peak_heights']
    peaks = peaks_obj[0]

    df = DataFrame({'peaks': peaks, 'heights': heights})
    idxmax = df['heights'].idxmax()

    idxs_remove = list(range(idxmax + 1))

    df = df.drop(idxs_remove)

    peaks_adj = df.nlargest(max_peak_count, ['heights'])

    return peaks_adj['peaks'].sort_values().to_numpy()


def generate_graph(arr: np.ndarray, max_diff: int) -> nx.DiGraph:
    G = nx.DiGraph()

    for p in arr:
        G.add_node(p)

    i = 0
    while i < arr.size:
        j = i + 1
        while j < arr.size:
            diff = arr[j] - arr[i]
            if diff <= max_diff:
                G.add_edge(arr[i], arr[j], length=diff)
            j += 1
        i += 1

    return G


def generate_combinations(G: nx.DiGraph, ladder_peak_count: int):

    start_nodes = [node for node in G.nodes if G.in_degree(node) == 0]
    end_nodes = [node for node in G.nodes if G.out_degree(node) == 0]

    if len(start_nodes) > 1:
        raise Exception("Can't generate. Too many start nodes.")

    if len(end_nodes) > 1:
        raise Exception("Can't generate, Too many end nodes.")

    all_paths = list(nx.all_simple_paths(G, start_nodes[0], end_nodes[0]))

    for p_arr in all_paths:
        for i in range(0, len(p_arr) - ladder_peak_count + 1):
            yield p_arr[i:i+ladder_peak_count]


# Below: Closure, class and function (for use with funtools.partial) for parallel processing of statistics.
#        All were slower than using a single thread.
#
# def make_statistics_fun(ladder):
#     def inner(combination):
#         c_arr = np.array(combination)
#
#         corr_peaks = stats.pearsonr(ladder, c_arr)
#
#         l_diff = np.diff(ladder)
#         c_diff = np.diff(c_arr)
#
#         corr_diffs = stats.pearsonr(l_diff, c_diff)
#
#         obj = {
#             'corr_peaks': corr_peaks.statistic,
#             'corr_diffs': corr_diffs.statistic,
#             'ladder': ladder.tolist(),
#             'peaks': c_arr.tolist()
#         }
#
#         return obj
#
#     return inner
#
#
# def get_statistics(ladder, combination):
#     c_arr = np.array(combination)
#
#     corr_peaks = stats.pearsonr(ladder, c_arr)
#
#     l_diff = np.diff(ladder)
#     c_diff = np.diff(c_arr)
#
#     corr_diffs = stats.pearsonr(l_diff, c_diff)
#
#     obj = {
#         'corr_peaks': corr_peaks.statistic,
#         'corr_diffs': corr_diffs.statistic,
#         'ladder': ladder.tolist(),
#         'peaks': c_arr.tolist()
#     }
#
#     return obj
#
#
# class Statistics:
#     def __init__(self, ladder):
#         self.ladder = ladder
#
#     def process(self, combination):
#         c_arr = np.array(combination)
#
#         corr_peaks = stats.pearsonr(self.ladder, c_arr)
#
#         l_diff = np.diff(self.ladder)
#         c_diff = np.diff(c_arr)
#
#         corr_diffs = stats.pearsonr(l_diff, c_diff)
#
#         obj = {
#             'corr_peaks': corr_peaks.statistic,
#             'corr_diffs': corr_diffs.statistic,
#             'ladder': self.ladder.tolist(),
#             'peaks': c_arr.tolist()
#         }
#
#         return obj


def main():
    # Testdata:
    # peaks = np.array([
    #     1577, 1735, 1845, 1895, 2060, 2222, 2386, 2551, 2668, 2717,
    #     2886, 2969, 3054, 3224, 3396, 3517, 3566, 3738, 3909, 4081, 4250, 4369, 4418, 4586,
    #     4750, 4913, 5073, 5182, 5226, 5380, 5528, 5671, 5810
    #     ])
    #

    ladder = np.array([
        100,
        114, 120, 140, 160, 180, 200,
        214, 220, 240, 260, 280, 300,
        314, 320, 340, 360, 380, 400,
        414, 420, 440, 460, 480, 500,
        514, 520, 540, 560, 580
    ])

    files = Path('demo/4071_Dx 230113_PRT1_PRT3_rn/').glob('*.fsa')

    # files = Path('demo/4062_Dx/').glob('*.fsa')
    raw_data = dict()

    for f in files:
        fn = f.name
        raw_data[fn] = np.array(SeqIO.read(f, 'abi').annotations["abif_raw"]['DATA205'])

    for sample in raw_data:
        print(f"processing {sample} ...")
        peaks = get_peaks(raw_data[sample], 38)

        diff = np.diff(peaks)
        max_diff = np.max(diff)

        max_len = max_diff * 1.5  # max length of peak jumps

        G = generate_graph(peaks, max_len)
        combs = list(generate_combinations(G, ladder.size))

        result = []

        for c in combs:
            c_arr = np.array(c)

            corr_peaks = stats.pearsonr(ladder, c_arr)

            l_diff = np.diff(ladder)
            c_diff = np.diff(c_arr)

            corr_diffs = stats.pearsonr(l_diff, c_diff)

            obj = {
                'corr_peaks': corr_peaks.statistic,
                'corr_diffs': corr_diffs.statistic,
                'ladder': ladder.tolist(),
                'peaks': c_arr.tolist()
            }

            result.append(obj)

        df = pd.DataFrame.from_records(result)
        df = df.sort_values(by='corr_peaks', ascending=False)
        best = df.iloc[0]

        print(f"results for {sample}")
        print(best)


if __name__ == '__main__':
    main()
