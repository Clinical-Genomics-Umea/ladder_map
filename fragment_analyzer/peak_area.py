import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.integrate import simpson
from lmfit.models import VoigtModel, GaussianModel
import lmfit
from collections import namedtuple
from typing import NamedTuple


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
