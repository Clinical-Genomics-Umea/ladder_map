import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
from lmfit.models import VoigtModel, GaussianModel, LorentzianModel
from fragment_analyzer.ladder_map import LadderMap


class PeakArea:
    def __init__(
        self,
        laddermap: LadderMap,
        model: str,
        channel: str = "DATA1",
        min_ratio: float = 0.2,
    ) -> None:
        self.file_name = laddermap.data_.parts[-1]
        self.raw_data = laddermap.adjusted_step_dataframe(channel=channel)

        # find peaks
        self.find_peaks_agnostic(min_ratio=min_ratio)

        # if no peaks could be found
        self.found_peaks = True
        if self.peak_information.shape[0] == 0:
            self.found_peaks = False
            print(
                f"No peaks could be found in {self.file_name}. Please look at the raw data."
            )

        # if peaks are found
        if self.found_peaks:
            print(f"{self.peak_information.shape[0]} peaks found in {self.file_name}")
            # find peak widths
            self.find_peak_widths()
            # divide peaks into individual dataframes
            self.divide_peaks()
            self.fit_df, self.fit_params, self.fit_report = self.fit_lmfit_model(
                model_=model
            )
            # calculate quotient
            self.calculate_quotient()

    def find_peak_widths(self, rel_height: float = 0.95):
        widths = peak_widths(
            self.peaks_dataframe.peaks,
            self.peaks_index,
            rel_height=rel_height,
        )

        df = pd.DataFrame(widths).T
        df.columns = ["x", "peak_height", "peak_start", "peak_end"]

        self.peak_widths = (
            df.assign(peak_start=lambda x: np.floor(x.peak_start).astype(int))
            .assign(peak_end=lambda x: np.ceil(x.peak_end).astype(int))
            .assign(peak_name=lambda x: range(1, x.shape[0] + 1))
        )

    @property
    def plot_peak_widths(self):
        fig = plt.figure(figsize=(20, 10))

        df = self.peaks_dataframe.loc[
            lambda x: x.step_adjusted > self.peak_information.step_adjusted.min() - 10
        ].loc[
            lambda x: x.step_adjusted < self.peak_information.step_adjusted.max() + 10
        ]

        plt.plot(df.step_adjusted, df.peaks)
        plt.plot(self.peak_information.step_adjusted, self.peak_information.peaks, "o")

        plt.ylabel("intensity")
        plt.xlabel("basepairs")
        plt.grid()

        return fig

    def divide_peaks(self, padding: int = 4):
        # add some padding to the left and right to be sure to include everything in the peak
        self.divided_peaks = [
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

        # if there only is 1 peak, return 1
        if len(areas) == 1:
            self.quotient = 1.0
            return

        # if there only are 2 peaks, return the quotient
        if len(areas) == 2:
            self.quotient = areas[1] / areas[0]
            return

        # return the last peak divided by the mean of the peaks to the left of it
        self.quotient = areas[-1] / areas[:-1].mean()

    @property
    def plot_lmfit_model(self):

        fig, axs = plt.subplots(1, len(self.fit_df), sharey=True, figsize=(20, 10))

        # if there is only one peak
        if len(self.fit_df) == 1:
            axs.plot(self.fit_df[0].step_adjusted, self.fit_df[0].peaks, "o")
            axs.plot(self.fit_df[0].step_adjusted, self.fit_df[0].fitted)
            axs.set_title(f"Peak 1 area: {self.fit_params[0]['amplitude']: .1f}")
            axs.grid()
        # if more than one peak
        else:
            for i, ax in enumerate(axs):
                ax.plot(self.fit_df[i].step_adjusted, self.fit_df[i].peaks, "o")
                ax.plot(self.fit_df[i].step_adjusted, self.fit_df[i].fitted)
                ax.set_title(
                    f"Peak {i + 1} area: {self.fit_params[i]['amplitude']: .1f}"
                )
                ax.grid()

        fig.suptitle(f"Quotient: {self.quotient: .2f}")
        fig.legend(["Raw data", "Model"])
        fig.supxlabel("basepairs")
        fig.supylabel("intensity")

        return fig

    @property
    def peak_position_area_dataframe(self) -> pd.DataFrame:
        """
        Returns a DataFrame of each peak and its properties
        """
        dataframes = []
        for i, _ in enumerate(self.fit_df):
            df = (
                self.fit_df[i]
                .loc[lambda x: x.peaks == x.peaks.max()]
                .assign(area=self.fit_params[i]["amplitude"])
                .assign(peak_name=f"Peak {i + 1}")
                .drop(columns="step_raw")
                .reset_index(drop=True)
                .rename(
                    columns={
                        "peaks": "peak_height",
                        "step_adjusted": "basepairs",
                        "fitted": "fitted_peak_height",
                    }
                )
                .drop_duplicates("peak_name")
                .assign(file_name=self.file_name)
            )
            dataframes.append(df)
        return pd.concat(dataframes).assign(quotient=self.quotient)

    # TODO
    # change peak_height to something appropriate... but what?
    # change min_ratio to something appropriate... but what?
    def find_peaks_agnostic(
        self, peak_height: int = 500, min_ratio: float = 0.2
    ) -> None:
        peaks_dataframe = self.raw_data.loc[lambda x: x.step_adjusted > 50]
        peaks_index, _ = find_peaks(peaks_dataframe.peaks, height=peak_height)

        peak_information = (
            peaks_dataframe.iloc[peaks_index]
            .assign(peaks_index=peaks_index)
            .assign(ratio=lambda x: x.peaks / x.peaks.max())
            .loc[lambda x: x.ratio > min_ratio]
            .assign(peak_name=lambda x: range(1, x.shape[0] + 1))
        )

        # update peaks_index based on the above filtering
        peaks_index = peak_information.peaks_index.to_numpy()

        # update class attributes
        self.peaks_index = peaks_index
        self.peaks_dataframe = peaks_dataframe
        self.peak_information = peak_information

    @property
    def plot_raw_data(self):
        """
        Plot the whole area of the raw data
        """
        fig = plt.figure(figsize=(20, 10))
        plt.plot(self.raw_data.step_adjusted, self.raw_data.peaks)

        plt.ylabel("intensity")
        plt.xlabel("basepairs")
        plt.grid()
        plt.title("Raw data")

        return fig
