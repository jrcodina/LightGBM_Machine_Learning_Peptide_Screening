# LightGBM Machine Learning Framework for Screening Large Peptide Libraries
Here we present a comprehensive guide on how to use Light Gradient Boost Machine, a machine learning framework, to screen large peptide libraries. By starting with a molecular docking of 1% of the peptide library, we are able to predict wich peptides are the best 5%. This reduces the time for docking in x20 by discarding 94% of the peptide library.  

## Pre-requesites
To follow this guide and use this tool you will need [R Project for Statistical Computing](https://www.r-project.org/) (version 4.2.2 (2022-10-31 ucrt) or higher). The R packages: [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html); [lightgbm](https://cran.microsoft.com/snapshot/2022-04-06/web/packages/lightgbm/index.html); [caret](https://cran.r-project.org/web/packages/caret/); [Peptides](https://cran.microsoft.com/snapshot/2022-01-28/web/packages/Peptides/index.html)

## First steps
The first thing you will need is to do the docking to 1% of the peptide library. Here we won't go into detail on how to perform moleculr docking. We have choosen the tetramer library, this is, a library containing all possible tetramers. If we consider only the 20 natural amino acids, then this library is 20^4 = 160,000 peptides long. 1% of this library is 1,600. We need to select 1,600 tetramers randomly. Then do molecular docking of these to a target. We can use a fast and rigid docking since we are doing already a large number of peptides. Our goal here is not to get an accurate binding score but rather to get a sense of what peptides will be better candidates than others. As we mentioned, this guide will not cover molecular docking. We will assume that step and start with that information. 

