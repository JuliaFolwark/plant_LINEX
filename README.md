# pLINEX: Lipid Network Explorer for Plants

## Abstract
This repository contains the implementation of the Lipid Network Explorer for plants (pLINEX), a bioinformatics tool developed to analyze stress responses in *Arabidopsis thaliana* by extending the existing Lipid Network Explorer (LINEX), originally designed for animals. Lipids are critical in plant metabolism, including maintaining cell structure, signaling, and responding to environmental stress. As climate change presents increasing challenges to agriculture, understanding plant lipidomes is essential for improving crop resilience. The primary objective of this research was to develop pLINEX, a plant-specific extension of LINEX, and to test it using lipidomics data from *A. thaliana* under stress conditions. While pLINEX successfully constructed plant-specific lipid networks, it did not identify significant changes in the lipidome under the various stress conditions, raising questions about the tool's effectiveness or the quality of the input data.

## Project Overview
pLINEX is a bioinformatics tool designed to analyze lipidomic data specific to plants. It extends the functionality of the original LINEX tool to accommodate the unique lipid structures and metabolic pathways found in plants. This tool was developed with a focus on *Arabidopsis thaliana*, a model organism in plant biology, and was tested using lipidomics data from experiments involving stress conditions.

### Key Features
- **Plant-Specific Lipid Network Construction**: pLINEX parses and curates data from the Plant Metabolic Network (PMN) and Rhea databases to build a plant-specific lipid network.
- **Stress Response Analysis**: The tool was developed to analyze changes in lipid composition and metabolism in various plant species under different stress conditions.
- **Customizable Data Input**: Users can input their lipidomics data, choose specific plant species, and select relevant databases for network construction.


### Usage

To utilize the content of this repository effectively, please follow these instructions:

**PMN Parser**:

* **Preperation:** Ensure you have the contents of the "PMN_parser" folder. Install the necessary dependencies listed in the "requirements.txt" file.
* **Execution**:Use the "pmn_file_parser.py" script to generate curated PMN database files. Please note that you will need to adjust the file paths within the script to match your local setup and provide the path to the database files. For guidance, refer to the annotations within the code.


**pLINEX**:

* **Importing adjusted pLINEX**: The package named "linex2" must be imported into your Python site-packages. This package includes the PMN database files. The PMN parser does not need to be executed beforehand unless you want to use a new version of the database, in which case you should update the package accordingly.
* **Using the pLINEX Web Tool**:Navigate to the linex package, which contains the modified web tool for pLINEX. Enter the "lipid_network_project" directory and run "start_linex.sh". Open your web browser and navigate to localhost:7000 to access the web tool. The web tool can then be used as described in the LINEX2 documentation. Note that only LINEX2 has been modified for pLINEX; the example data has not been updated and should not be used. If you need to close the web application, ensure that all processes associated with LINEX are terminated. Run stop_linex.sh to identify and stop the relevant processes. This step is crucial to ensure that the tool can be restarted without issues in future sessions.


### Publication
If you use this package, please cite, since this package is based on the LINEX2 tool:

Rose and Koehler et al. "**Lipid network and moiety analysis for revealing enzymatic dysregulation and mechanistic alterations from lipidomics data**",
_Briefings in Bioinformatics_ **2023**, bbac572; doi: [https://doi.org/10.1093/bib/bbac572](https://doi.org/10.1093/bib/bbac572)

### Data in package

The package includes data from the [Reactome](https://reactome.org/) and 
[Rhea](https://www.rhea-db.org/) databases.

* **Rhea** data is available under the 
[Creative Commons Attribution (CC BY 4.0) License](https://creativecommons.org/licenses/by/4.0/)
and is published here: Bansal et al. "Rhea, the reaction knowledgebase in 2022".
Nucleic Acids Res. (2021), DOI: [10.1093/nar/gkab1016](https://doi.org/10.1093/nar/gkab1016)
* **Reactome** data is available under the
[Creative Commons Public Domain (CC0) License](https://creativecommons.org/publicdomain/zero/1.0/)
and is published here: Jassal et al. "The reactome pathway knowledgebase".
Nucleic Acids Res. (2020), DOI: [10.1093/nar/gkz1031](https://doi.org/10.1093/nar/gkz1031)
* **Plant Metabolic Network** data is available under the
[Creative Commons Public Domain (CC0) License](https://creativecommons.org/publicdomain/zero/1.0/)
and is published here: Hawkins et al. "Plant Metabolic Network 15: A resource of genome-wide metabolism databases for 126 plants and algae".
Journal of Integrative Plant Biology (2021), DOI: [10.1111/jipb.13163](https://doi.org/10.1111/jipb.13163)


### License

The LINEX tool was published under the AGPLv3 license.

![AGPLv3 logo](https://www.gnu.org/graphics/agplv3-with-text-162x68.png)

The software includes code from the 
* **Rhea** data is available under the 
[Creative Commons Attribution (CC BY 4.0) License](https://creativecommons.org/licenses/by/4.0/)
and is published here: Bansal et al. "Rhea, the reaction knowledgebase in 2022".
Nucleic Acids Res. (2021), DOI: [10.1093/nar/gkab1016](https://doi.org/10.1093/nar/gkab1016)
* **Reactome** data is available under the
[Creative Commons Public Domain (CC0) License](https://creativecommons.org/publicdomain/zero/1.0/)
and is published here: Jassal et al. "The reactome pathway knowledgebase".
Nucleic Acids Res. (2020), DOI: [10.1093/nar/gkz1031](https://doi.org/10.1093/nar/gkz1031)
* **Plant Metabolic Network** data is available under the
[Open Database (ODbL) License](https://opendatacommons.org/licenses/odbl/)
and is published here: Hawkins et al. "Plant Metabolic Network 15: A resource of genome-wide metabolism databases for 126 plants and algae".
Journal of Integrative Plant Biology (2021), DOI: [10.1111/jipb.13163](https://doi.org/10.1111/jipb.13163)


