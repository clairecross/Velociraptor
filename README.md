# Velociraptor

### Velociraptor publication available through bioRxiv: https://www.biorxiv.org/content/10.1101/2024.05.01.591375v1

Velociraptor seeks cells based on a user-defined phenotype. With VR-Eye, cell identity can be defined based on phenotypic comparisons with well-established cell identities to assess if a population strongly matches a single cell type or if it has similarity with multiple populations. VR-Eye also enables automated cross-platform cell identification.

### Figure 1:

![alt text](https://github.com/clairecross/Velociraptor/blob/main/Velociraptor%20Overview.png)

Velociraptor Overview. a) Velociraptor seeks cells based on a user-defined reference phenotype using MEM label syntax where “+0” indicates a feature that is specifically lacking from a target population (e.g., no expression in that population) and “+10” indicates a feature that is specifically present in that population (e.g., high, uniform expression in all cells of that population). The phenotype for every cell in a dataset is calculated using that cell and its closest phenotypic neighbors (CELL), identified using k-Nearest Neighbors (KNN).  The difference between the sought reference phenotype (REF) and the cell’s phenotype is calculated using root mean square deviation (RMSD) and recorded as a per-cell feature. This is repeated for each reference population supplied as input. b) Velociraptor similarity scores can be visualized per-cell and plotted in the coordinate space of the original image. In this example, cells are shaded according to their degree of matching the specified B cell reference population. A spectrum intensity scale indicates cell similarity with red indicating high similarity and purple indicating low similarity.

Velociraptor was developed in the laboratory of Dr. Jonathan Irish at Vanderbilt University. The research was supported by the following funding resources: NIH/NCI grants R01 NS096238 (RAI, JMI), R01 CA226833 (JMI, CEC, SM, MJH), R01 NS118580 (RAI), U01 AI125056 (JMI), U54 CA217450 (JMI, MJH), T32GM137793 (CEC), the Vanderbilt-Ingram Cancer Center (VICC, P30 CA68485), the Michael David Greene Brain Cancer Fund (RAI, JMI), the Southeastern Brain Tumor Foundation (RAI, JMI), a gift from Daniel F Hewins (RAI), the Ben & Catherine Ivy Foundation (RAI, JMI), and by the Human Immunology Discovery Initiative of the Vanderbilt Center for Immunobiology.

If you’re interested in learning more, check out the other tools on the CytoLab Github page at:
https://github.com/cytolab/

### Datasets used here:
1. The human glioblastoma mass cytometry dataset (RAPID dataset) was downloaded from: http://flowrepository.org/id/FR-FCM-Z24K. This study can be found at: Leelatian, N., Sinnaeve, J. et al. Unsupervised machine learning reveals risk stratifying glioblastoma tumor cells. eLife 9 (2020). https://doi.org/10.7554/eLife.56879.

2. The human breast cancer imaging mass cytometry dataset (IMC dataset) was downloaded from: Tietscher Sandra. (2022). Imaging Mass Cytometry Dataset of exhausted and non-exhausted breast cancer microenvironments [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4911135. This study can be found at: Tietscher, S., Wagner, J., Anzeneder, T. et al. A comprehensive single-cell map of T cell exhaustion-associated immune environments in human breast cancer. Nat Commun 14, 98 (2023). https://doi.org/10.1038/s41467-022-35238-w.

