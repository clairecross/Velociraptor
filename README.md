# Velociraptor

### Velociraptor publication available through bioRxiv: https://www.biorxiv.org/content/10.1101/2024.05.01.591375v1

The Velociraptor workflow includes two novel cell identification tools: Velociraptor-Claw (VR-Claw) and Velociraptor-Eye (VR-Eye).

VR-Claw reveals condition-specific cell populations (e.g., cells associated with overall survival).

VR-Eye seeks cells based on a user-defined phenotype. With VR-Eye, cell identity can be defined based on phenotypic comparisons with well-established cell identities to assess if a population strongly matches a single cell type or if it has similarity with multiple populations. VR-Eye also enables automated cross-platform cell identification.

### Figure 1:

![alt text](https://github.com/clairecross/Velociraptor/blob/main/Velociraptor%20Overview.png)

Velociraptor Overview. VR-Claw utilizes local cell neighborhoods to reveal clinically relevant cell populations. Local cell neighborhoods are defined by KNN (depicted in light grey) and tested for association with patient survival. Prognostic populations are clustered and described using Marker Enrichment Modeling (MEM). MEM labels are used to quantitatively define reference phenotypes that are sought in new samples with VR-Eye. Each cell’s similarity to the reference phenotype is calculated and visualized where red indicates high similarity and blue indicates low similarity. Multiple phenotypes can be sought and used to build a multidimensional model of cell identity. 

Velociraptor was developed in the laboratory of Dr. Jonathan Irish at Vanderbilt University.  The research was supported by the following funding resources: NIH/NCI grants R01 NS096238 (RAI, JMI), R01 CA226833 (JMI, CEC, SM, MJH), R01 NS118580 (RAI), U01 AI125056 (JMI), U54 CA217450 (JMI, MJH), T32GM137793 (CEC), the Vanderbilt-Ingram Cancer Center (VICC, P30 CA68485), the Michael David Greene Brain Cancer Fund (RAI, JMI), the Southeastern Brain Tumor Foundation (RAI, JMI), a gift from Daniel F Hewins (RAI), the Ben & Catherine Ivy Foundation (RAI, JMI), and by the Human Immunology Discovery Initiative of the Vanderbilt Center for Immunobiology.

If you’re interested in learning more, check out the other tools on the CytoLab Github page at:
https://github.com/cytolab/

### Datasets used here:
1. The human glioblastoma mass cytometry dataset (RAPID dataset) was downloaded from: http://flowrepository.org/id/FR-FCM-Z24K. This study can be found at: Leelatian, N., Sinnaeve, J. et al. Unsupervised machine learning reveals risk stratifying glioblastoma tumor cells. eLife 9 (2020). https://doi.org/10.7554/eLife.56879.

2. The human breast cancer imaging mass cytometry dataset (IMC dataset) was downloaded from: Tietscher Sandra. (2022). Imaging Mass Cytometry Dataset of exhausted and non-exhausted breast cancer microenvironments [Data set]. Zenodo. https://doi.org/10.5281/zenodo.4911135. This study can be found at: Tietscher, S., Wagner, J., Anzeneder, T. et al. A comprehensive single-cell map of T cell exhaustion-associated immune environments in human breast cancer. Nat Commun 14, 98 (2023). https://doi.org/10.1038/s41467-022-35238-w.

