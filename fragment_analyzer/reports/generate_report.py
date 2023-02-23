import panel as pn
from pathlib import Path

import fragment_analyzer

pn.extension("tabulator")
pn.extension("vega", sizing_mode="stretch_width", template="fast")
pn.widgets.Tabulator.theme = "modern"


def generate_report(report_type: str, fsa_file: str, ladder: str, folder: str) -> None:
    """
    Generates the report
    """
    match report_type:
        case "peak_area_report":
            peak_area_report(fsa_file, ladder, folder)
        case _:
            print("No support")



def header(
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


def generate_peak_area_report(name: str, plot_raw, plot_ladder, plot_peaks, peak_area):
    head = header(
        text=f"""
        # Fragment Analysis Report
        ## Report of {name}
        """,
        fontsize="20px",
        bg_color="#03a1fc",
        height=185,
    )
    ### ----- Raw Data plot ----- ###
    raw_data_markdown = header("# Raw Data plot", height=100)
    raw_data_plot = pn.pane.Matplotlib(plot_raw.plot_raw_data)

    ### ----- Ladder info and raw data----- ###
    best_ladder_markdown = header("# Fit of the Ladder", height=100)
    best_ladder_plot = pn.pane.Matplotlib(plot_ladder.plot_ladder_peaks)

    correlation_plot = pn.pane.Matplotlib(
        plot_ladder.plot_model_fit,
    )

    ### ----- Peaks info ----- ###
    peaks_markdown = header("# Peaks", height=100)
    peaks_plot = pn.pane.Matplotlib(plot_peaks.plot_peaks)

    ### ----- Quotient info ----- ###
    quotient_markdown = header("# Areas", height=100)
    quotient_plot = pn.pane.Matplotlib(plot_peaks.plot_areas)

    ### ----- Quotient info ----- ###
    df_markdown = header("# Peak Information Table", height=100)
    df_table = pn.widgets.Tabulator(
        peak_area.peak_position_area_dataframe,
        show_index=False,
        name="Peak information",
    )

    ### ----- Model fitting info ----- ###
    model_fitting_markdown = header("# Fitting of the Models", height=100)
    model_fitting_report = []
    for i, x in enumerate(peak_area.fit_report):
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


def generate_peak_area_no_peaks(name, plot_raw):
    head = header(
        text=f"""
        # Fragment Analysis Report
        ## Report of {name}
        """,
        fontsize="20px",
        bg_color="#03a1fc",
        height=185,
    )

    no_peaks_markdown = header(
        "# No peaks could be generated. Please look at the raw data.", height=100
    )
    raw_plot = pn.pane.Matplotlib(plot_raw.plot_raw_data)
    return pn.Column(
        head,
        no_peaks_markdown,
        raw_plot,
    )



def peak_area_report(fsa_file: str,ladder: str, folder: str) -> None:
    """
    """
    fsa = fragment_analyzer.FsaFile(fsa_file, ladder)
    file_name = fsa.file_name
    ladder_assigner = fragment_analyzer.PeakLadderAssigner(fsa)
    model = fragment_analyzer.FitLadderModel(fsa, ladder_assigner)
    raw_plots = fragment_analyzer.PlotRawData(model)
    ladder_plots = fragment_analyzer.PlotLadder(model)
    peak_areas = fragment_analyzer.PeakArea(model, peak_finding_model="gauss")
    peak_plots = fragment_analyzer.PlotPeakArea(peak_areas)

    # create the output folder if it doesn't exist
    outpath = Path(folder)
    if not outpath.exists():
        outpath.mkdir(parents=True)

    # If no peaks could be found
    if not peak_areas.found_peaks:
        outname = outpath / f"FAILED-fragment_analysis-report-{file_name}.html"
        generate_peak_area_no_peaks(file_name, raw_plots).save(
            outname, title=file_name
        )

    else:
        outname = outpath / f"fragment_analysis-report-{file_name}.html"
        generate_peak_area_report(
            file_name, raw_plots, ladder_plots, peak_plots, peak_areas
        ).save(outname, title=file_name)
