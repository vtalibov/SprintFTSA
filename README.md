# SprintFTSA

Script for bulk processing of thermal shift assay (FTSA/DSF) data. Uses a simple GUI to prompt for input and output files.

![](/SprintFTSA.png)

### Data

Written for .csv export of data from Biorad RT-PCR thermo cyclers, *e.g.* CFX96/CFX384; see input data example.

### Data analysis

Transition phase of the thermal denaturation curve and the corresponding temperature range are identified numerically, then the curve is approximated with a 5-parametric logistic curve equation. The script outputs results table and the corresponding plots.

![](/results_example/C4.png)

### Good to know

Samuel, E. L. G.; Holmes, S. L.; Young, D. W. [Processing Binding Data Using an Open-Source Workflow.](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00577-1) J. Cheminform. 2021, 13 (1), 99. **Good reference on non-mechanistic curve fitting of FTSA data.**