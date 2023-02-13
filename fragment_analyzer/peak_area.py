import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
from lmfit.models import VoigtModel, GaussianModel, LorentzianModel
import lmfit
from collections import namedtuple
from typing import NamedTuple


class PeakArea:
    def __init__(
        self,
        data: pd.DataFrame,
        start: int,
        end: int,
        num_peaks: int,
        model: str,
        rel_height: float = 0.95,
        padding: int = 3,
    ) -> None:
        self.data = data
        self.peaks_index, self.peaks_dataframe = self.find_fragment_peaks(
            start=start, end=end, num_peaks=num_peaks
        )
        self.peak_widths = self.find_fragment_peak_width(rel_height=rel_height)
        self.divided_peaks = self.divide_peaks(padding=padding)
        self.fit_df, self.fit_params, self.fit_report = self.fit_lmfit_model(
            model_=model
        )
        self.quotient = self.calculate_quotient()

    def find_fragment_peaks(self, start, end, num_peaks: int, threshold: int = 10000):
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
                    f"No peaks found in {start=}, {end=}. Try another part of the graph."
                )

        return peaks, data

    def find_fragment_peak_width(self, rel_height: float = 0.95):
        widths = peak_widths(
            self.peaks_dataframe.peaks,
            self.peaks_index,
            rel_height=rel_height,
        )

        df = pd.DataFrame(widths).T
        df.columns = ["x", "peak_height", "peak_start", "peak_end"]

        return (
            df.assign(peak_start=lambda x: np.floor(x.peak_start).astype(int))
            .assign(peak_end=lambda x: np.ceil(x.peak_end).astype(int))
            .assign(peak_name=lambda x: range(1, x.shape[0] + 1))
        )

    def plot_peak_widths(self):
        fig = plt.figure(figsize=(15, 10))

        plt.plot(self.peaks_dataframe.step_adjusted, self.peaks_dataframe.peaks)

        plt.plot(
            self.peaks_dataframe.step_adjusted.iloc[self.peaks_index],
            self.peaks_dataframe.peaks.iloc[self.peaks_index],
            "o",
        )

        for widths in self.peak_widths.itertuples():
            plt.hlines(
                y=widths.peak_height,
                xmin=self.peaks_dataframe.step_adjusted.iloc[widths.peak_start],
                xmax=self.peaks_dataframe.step_adjusted.iloc[widths.peak_end],
                color="C3",
            )
        plt.ylabel("intensity")
        plt.xlabel("basepairs")
        plt.grid()

        return fig

    def divide_peaks(self, padding):
        # add some padding (defalt 3) to the left and right to be sure to include everything in the peak
        return [
            self.peaks_dataframe.iloc[x.peak_start - padding : x.peak_end + padding]
            for x in self.peak_widths.itertuples()
        ]

    def fit_lmfit_model(self, model_: str):
        if model_ == "gauss":
            model = GaussianModel()
        elif model_ == "voigt":
            model = VoigtModel()
        elif model_ == "lorentzian":
            model = LorentzianModel()
        else:
            raise NotImplementedError(
                f"{model_} is not implemented! Options: [gauss, voigt, lorentzian]"
            )

        fitted_df = []
        fitted_parameters = []
        fitted_report = []
        for df in self.divided_peaks:
            df = df.copy()
            y = df.peaks.to_numpy()
            x = df.step_adjusted.to_numpy()

            params = model.guess(y, x)
            out = model.fit(y, params, x=x)

            fitted_df.append(df.assign(fitted=out.best_fit, model=model_))
            fitted_parameters.append(out.values)
            fitted_report.append(out.fit_report())

        return fitted_df, fitted_parameters, fitted_report

    def calculate_quotient(self):
        areas = np.array([x["amplitude"] for x in self.fit_params])

        # if there only are 2 peaks, return the quotient
        if len(areas) == 2:
            return areas[1] / areas[0]

        # return the last peak divided by the mean of the peaks to the left of it
        return areas[-1] / areas[:-1].mean()

    def plot_lmfit_model(self):

        fig, axs = plt.subplots(1, len(self.fit_df), sharey=True, figsize=(15, 10))

        for i, ax in enumerate(axs):
            ax.plot(self.fit_df[i].step_adjusted, self.fit_df[i].peaks, "o")
            ax.plot(self.fit_df[i].step_adjusted, self.fit_df[i].fitted)
            ax.set_title(f"Peak {i + 1} area: {self.fit_params[i]['amplitude']: .1f}")
            ax.grid()

        fig.suptitle(f"Quotient: {self.quotient: .2f}")
        fig.legend(["Raw data", "Model"])
        fig.supxlabel("basepairs")
        fig.supylabel("intensity")

        return fig
