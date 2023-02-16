import panel as pn
from pathlib import Path
from fragment_analyzer.ladder_map import LadderMap
from fragment_analyzer.peak_area import PeakArea

pn.extension("tabulator")
pn.extension("vega", sizing_mode="stretch_width", template="fast")
pn.widgets.Tabulator.theme = "modern"


class Report:
    def __init__(self, laddermap: LadderMap, peakarea: PeakArea):
        self.laddermap = laddermap
        self.peakarea = peakarea
        self.name = self.peakarea.file_name

    def header(
        self,
        text: str,
        bg_color: str = "#04c273",
        height: int = 150,
        fontsize: str = "px20",
        textalign: str = "center",
    ):
        """
        Template for markdown header like block
        """
        return pn.pane.Markdown(
            f"""
            {text}
            """,
            background=bg_color,
            height=height,
            margin=10,
            style={
                "color": "white",
                "padding": "10px",
                "text-align": f"{textalign}",
                "font-size": f"{fontsize}",
            },
        )

    def generate_report(self):
        head = self.header(
            text=f"""
            # Fragment Analysis Report
            ## Report of {self.name}
            """,
            fontsize="20px",
            bg_color="#03a1fc",
            height=185,
        )
        ### ----- Raw Data plot ----- ###
        raw_data_markdown = self.header("# Raw Data plot", height=100)
        raw_data_plot = pn.pane.Matplotlib(self.peakarea.plot_raw_data)

        ### ----- Ladder info and raw data----- ###
        best_ladder_markdown = self.header("# Fit of the Ladder", height=100)
        best_ladder_plot = pn.pane.Matplotlib(self.laddermap.plot_best_sample_ladder)

        correlation_plot = pn.pane.Matplotlib(
            self.laddermap.plot_ladder_correlation,
        )

        ### ----- Peaks info ----- ###
        peaks_markdown = self.header("# Peaks", height=100)
        peaks_plot = pn.pane.Matplotlib(self.peakarea.plot_peak_widths)

        ### ----- Quotient info ----- ###
        quotient_markdown = self.header("# Areas", height=100)
        quotient_plot = pn.pane.Matplotlib(self.peakarea.plot_lmfit_model)

        ### ----- Quotient info ----- ###
        df_markdown = self.header("# Peak Information Table", height=100)
        df_table = pn.widgets.Tabulator(
            self.peakarea.peak_position_area_dataframe,
            show_index=False,
            name="Peak information",
        )

        ### ----- Model fitting info ----- ###
        model_fitting_markdown = self.header("# Fitting of the Models", height=100)
        model_fitting_report = []
        for i, x in enumerate(self.peakarea.fit_report):
            markdown = pn.pane.Markdown(f"## Peak number {i + 1}:")
            model_fitting_report.append(markdown)
            model_fitting_report.append(x)
            model_fitting_report.append(pn.layout.Divider())

        ### CREATE REPORT ###

        return pn.Column(
            head,
            raw_data_markdown,
            raw_data_plot,
            pn.layout.Divider(),
            best_ladder_markdown,
            best_ladder_plot,
            correlation_plot,
            pn.layout.Divider(),
            peaks_markdown,
            peaks_plot,
            pn.layout.Divider(),
            quotient_markdown,
            quotient_plot,
            pn.layout.Divider(),
            df_markdown,
            df_table,
            pn.layout.Divider(),
            model_fitting_markdown,
            *model_fitting_report,
        )

    def generate_no_peaks_report(self):
        head = self.header(
            text=f"""
            # Fragment Analysis Report
            ## Report of {self.name}
            """,
            fontsize="20px",
            bg_color="#03a1fc",
            height=185,
        )

        no_peaks_markdown = self.header(
            "# No peaks could be generated. Please look at the raw data.", height=100
        )
        raw_plot = pn.pane.Matplotlib(self.peakarea.plot_raw_data)
        return pn.Column(
            head,
            no_peaks_markdown,
            raw_plot,
        )


def generate_report(laddermap: LadderMap, peakarea: PeakArea, folder: str) -> None:
    """
    Generates an HTML report for a given ladder map and peak area, and saves it to the specified folder.

    Args:
        laddermap: A LadderMap object containing the ladder map data.
        peakarea: A PeakArea object containing the peak area data.
        folder: A string representing the folder where the report will be saved.

    Returns:
        None

    Example usage:
    # create a LadderMap and PeakArea object
    laddermap = LadderMap(...)
    peakarea = PeakArea(...)

    # generate a report and save it to a folder called 'reports'
    generate_report(laddermap, peakarea, 'reports')
    """

    report = Report(laddermap, peakarea)

    # create the output folder if it doesn't exist
    outpath = Path(folder)
    if not outpath.exists():
        outpath.mkdir(parents=True)

    # If no peaks could be found
    if not peakarea.found_peaks:
        outname = outpath / f"FAILED-fragment_analysis-report-{report.name}.html"
        report.generate_no_peaks_report().save(outname, title=report.name)

    else:
        outname = outpath / f"fragment_analysis-report-{report.name}.html"
        report.generate_report().save(outname, title=report.name)
