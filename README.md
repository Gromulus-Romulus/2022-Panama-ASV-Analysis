# 2022-Panama-ASV-Analysis

The code here analyzes data from an ASV table of the ITS gene region.
Sequence variants and sample metadata were provided by the McGuire Lab at the University of Oregon (IE2).

`phyloseq` package was used to handle ASV data. `./Data/setup.R` imports all data and writes them to `.Rds` files, which are easily importable into other r scripts.

Other dependencies: `vegan`, `spaa`, and `EcolUtils` were used for ordination and PERMANOVA tests.

Each of the `.Rmd` files do a separate statistical analysis:
- `funguild.charts.Rmd` visualizes relative abundances of ecological guilds (determined prior via FunGUILD classification).
- `genera.charts.Rmd` visualizes relatie abundances of fungal genera.
- `ordination.plots.Rmd` uses non-metric multidimensional scaling (NMDS) to produce ordination plots of community divergence.
- `permanova.Rmd` uses a pairwise PERMANOVA test to confirm community divergence. Assumptions are tested using Beta-dispersion.

## Installing Dependencies (Updated March 16th, 2023):

Prior to use, you will need to have some packages pre-installed in order to run this code properly. Packages such as `vegan` and `spaa` may be downloaded in the usual way: using the `install.packages` function.

However, `phyloseq` and `EcolUtils` require special instructions for installation. You may find the links to these instructions below:

* Installing `phyloseq`: https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
* Installing `EcolUtils`: https://github.com/GuillemSalazar/EcolUtils

You may also want to consider installing RStudio, since most of the code files are written in R markdown (`.Rmd` format).
