# SprintFTSA

Script for bulk processing of thermal shift assay (FTSA/DSF) data. Uses a simple GUI to prompt for input and output files.

![](/SprintFTSA.png)

### Data

Written for .csv export of data from Biorad RT-PCR thermo cyclers, *e.g.* CFX96/CFX384; see input data example.

### Data analysis

Transition phase of the thermal denaturation curve and the corresponding temperature range are identified numerically, then the curve is approximated with a 5-parametric logistic curve equation. The script outputs results table and the corresponding plots.

![](/results_example/B2.png)

### Isothermal analysis

An optional feature, allows to roughly estimate unfolded fraction at specified temperature. Useful whenever FTSA have to be compared to quantitative CETSA or other thermal denaturation-based target engagement assays.

*Considering relatively high protein concentration in FTSA experiments, to quantify $$EC_{50}$$ analyze concentration-response curves with quadratic binding equations - tight-binding conditions are often satisfied.*

### Good to know

Samuel, E. L. G.; Holmes, S. L.; Young, D. W. [Processing Binding Data Using an Open-Source Workflow.](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00577-1) J. Cheminform. 2021, 13 (1), 99. **Good reference on non-mechanistic curve fitting of FTSA data.**