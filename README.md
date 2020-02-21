This file is a part of Di-PASS (DIP® Pathway Activation Scoring System)
The package is distributed under Di-PASS
Version 0.0.1

Copyright © 2020 Deep Intelligent Pharma.

Room 1608, Tower B T3, Wangjing SOHO, #1 Futong East Avenue,
Wangjing, Chaoyang District, Beijing, China

Di-PASS is a improved pathway analysis method that is based on iPANDA algorithm.
The method introduced here includes more elaborate KEGG pathway analysis and
gene importance estimation to robustly identify relevant pathways.

The introduction of the codes supplied within this method is described hereafter.

The codes are written in R and Python 3.6 The R code depends on KEGGgraph and RBGL
packages from Bioconductor, which can be obtained by the following links:

http://www.bioconductor.org/packages/release/bioc/html/KEGGgraph.html
http://www.bioconductor.org/packages/release/bioc/html/RBGL.html

The package was tested with the following configuration:

system: windows 7/10 professional
python version: 2.7.6
R version: 3.5.1
KEGGgraph version: 1.46.0
RBGL version: 1.62.1

1. pathway_topo.R
   Code is used to analysis pathway xml file. Files in file folder named "pathway_topo_input"  
   are xml files downloaded form KEGG database. The output files are in "process_topo_input" folder.

2. process_topo_kgml.py
   Code is used to calculate topology coefficients. Input file of the code contains three parts: output file
   from "pathway_topo.R", human gene symbol from NCBI, and coexpression units from cluster. The output files are input for iPANDA algorithm.

The details of the iPANDA algorithm can be found in the iPANDA manuscript:

    Ozerov I.V. et al., 2016, Nature Communications, 'In silico Pathway
    Activation Network Decomposition Analysis (iPANDA) as a method for biomarker
    development'.

Please cite the paper in any publication containing results acquired
using this method.