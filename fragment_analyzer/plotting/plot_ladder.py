import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from fragment_analyzer.ladder_fitting.fit_ladder_model import FitLadderModel


class PlotLadder:
    def __init__(self, model: FitLadderModel) -> None:
        self.model = model

    @property
    def plot_ladder_peaks(self) -> matplotlib.figure.Figure:
        trace = self.model.fsa_file.size_standard
        best_combination = self.model.best_combination
        ladder_name = self.model.fsa_file.ladder
        ladder_size = self.model.fsa_file.ref_sizes

        fig_ladder_peaks = plt.figure(figsize=(20, 10))
        plt.plot(trace)
        plt.plot(best_combination, trace[best_combination], "o")
        plt.xlabel("Time")
        plt.ylabel("Intensity")
        plt.legend(["Trace", "Peak (bp)"])
        plt.title(ladder_name)
        plt.grid()

        for peak, ladder in zip(best_combination, ladder_size):
            plt.text(peak, trace[peak], ladder)

        return fig_ladder_peaks

    @property
    def plot_model_fit(self):

        ladder_size = self.model.fsa_file.ref_sizes
        best_combination = self.model.best_combination

        predicted = self.model.model.predict(best_combination)
        ladder_name = self.model.fsa_file.ladder

        mse = self.model.mse
        r2 = self.model.r2

        fig_model_fit = plt.figure(figsize=(20, 10))
        plt.plot(ladder_size, best_combination, "o")
        plt.plot(predicted, best_combination, "x")
        plt.xticks(np.arange(0, np.max(ladder_size), 50))
        plt.xlabel("bp")
        plt.yticks(np.arange(0, np.max(best_combination), 500))
        plt.suptitle(ladder_name)
        plt.title(f"{mse=}, {r2=}")
        # plt.title(f"{self.model.model[0]=}")
        plt.legend(["True value", "Predicted value"])
        plt.grid()

        return fig_model_fit
