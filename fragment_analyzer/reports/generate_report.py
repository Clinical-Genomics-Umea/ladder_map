import pandas as pd
import numpy as np
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
        self.name = self.laddermap.data_.parts[-1]
        
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
        
        ### ----- Ladder info ----- ###
        best_ladder_markdown = self.header("# Fit of the Ladder", height=100)
        best_ladder_plot = pn.pane.Matplotlib(self.laddermap.plot_best_sample_ladder(), align="center")
    
        correlation_plot = pn.pane.Matplotlib(
            self.laddermap.plot_ladder_correlation(),
        )
        
        ### ----- Peaks info ----- ###
        peaks_markdown = self.header("# Peaks", height=100)
        peaks_plot = pn.pane.Matplotlib(
            self.peakarea.plot_peak_widths(), 
            align="center"
        )
        
        
        ### ----- Quotient info ----- ###
        quotient_markdown = self.header("# Areas", height=100)
        quotient_plot = pn.pane.Matplotlib(self.peakarea.plot_lmfit_model())
        
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
            model_fitting_markdown,
            *model_fitting_report
        )
    

def generate_report(laddermap: LadderMap, peakarea: PeakArea, folder: str):
    report = Report(laddermap, peakarea)
    
    outpath = Path(folder)
    if not outpath.exists():
        outpath.mkdir(parents=True)
        
    outname = outpath / f"fragment_analysis-report-{report.name}.html"
    report.generate_report().save(outname, title=report.name)
    