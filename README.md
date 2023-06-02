# LightGBM Machine Learning Framework for Screening Large Peptide Libraries
Here we present a comprehensive guide on how to use Light Gradient Boost Machine, a machine learning framework, to screen large peptide libraries. By starting with a molecular docking of 1% of the peptide library, we are able to predict wich peptides are the best 5%. This reduces the time for docking in x20 by discarding 94% of the peptide library.  

## Pre-requesites
To follow this guide and use this tool you will need [R Project for Statistical Computing](https://www.r-project.org/) (version 4.2.2 (2022-10-31 ucrt) or higher). Also the packages: [dplyr package](https://cran.r-project.org/web/packages/dplyr/index.html); [LightGBM package](https://cran.microsoft.com/snapshot/2022-04-06/web/packages/lightgbm/index.html); [caret package](https://cran.r-project.org/web/packages/caret/)
