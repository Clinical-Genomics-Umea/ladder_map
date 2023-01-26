## ladder map

Matches ladders to peaks by correlation for fragment analysis. The strategy resembles the one used by [Fragman](https://cran.r-project.org/web/packages/Fragman/index.html) for R.

One difference is that combinations of peaks are generated using [Networkx](https://networkx.org/) to eliminate impossible combinations. This reduces complexity substantially and allows for a greedy search to identify the best match.
