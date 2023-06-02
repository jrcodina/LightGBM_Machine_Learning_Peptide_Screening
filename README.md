# LightGBM Machine Learning Framework for Screening Large Peptide Libraries
Here we present a comprehensive guide on how to use Light Gradient Boost Machine, a machine learning framework, to screen large peptide libraries. By starting with a molecular docking of 1% of the peptide library, we are able to predict wich peptides are the best 5%. This reduces the time for docking in x20 by discarding 94% of the peptide library.  

## Pre-requesites
To follow this guide and use this tool you will need [R Project for Statistical Computing](https://www.r-project.org/) (version 4.2.2 (2022-10-31 ucrt) or higher). The R packages: [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html); [lightgbm](https://cran.microsoft.com/snapshot/2022-04-06/web/packages/lightgbm/index.html); [caret](https://cran.r-project.org/web/packages/caret/); [Peptides](https://cran.microsoft.com/snapshot/2022-01-28/web/packages/Peptides/index.html)

## First steps
The first thing you will need is to do the docking to 1% of the peptide library. Here we won't go into detail on how to perform moleculr docking. We have choosen the tetramer library, this is, a library containing all possible tetramers. If we consider only the 20 natural amino acids, then this library is 20^4 = 160,000 peptides long. 1% of this library is 1,600. We need to select 1,600 tetramers randomly. Then do molecular docking of these to a target. We can use a fast and rigid docking since we are doing already a large number of peptides. Our goal here is not to get an accurate binding score but rather to get a sense of what peptides will be better candidates than others. As we mentioned, this guide will not cover molecular docking. We will assume that step and start with that information. 

# (PART) Prepare Data {-}

First we will extract the features of the full tetramer library which is the library we want to screen. 

## Load libraries
```{r, eval=F, libraries1}

if (!require(dplyr)) {
  install.packages("dplyr")
}

if (!require(Peptides)) {
  install.packages("Peptides")
}


library(dplyr)
library(Peptides)
```

## Obtain properties
First we obtain the properties for the tetramer library. The function 'PCAprop' can be found at the end of this section. It extracts 99 properties form the [Peptides](https://cran.microsoft.com/snapshot/2022-01-28/web/packages/Peptides/index.html)  library for all tetramers. The function generates all posible combinaitons of the 20 natural amino acids for each position. If n=4 it will generate 20^4 = 160,000 peptides and their 99 properties. Take into account that this scalates exponentially, and a matrix of hexamers (64 Million peptide x 99 properties) can take along time to be generated and take all your ram. 

```{r,  eval=F, properties}
#Full library
tetra_160 <- PCAprop(n = 4)
```

## Normalize
Then we normalize the full tetra-peptide library. Eventhough LighGBM doesn't require normalization, this step will help us if we want to test other algorithms. 
```{r,  eval=F, normalize}
#Min-Max normalization
normalize <- function(x){
X <- x[,-1]
Y <- x[,1]
  
max = apply(X , 2 , max)
min = apply(X, 2 , min)
Xn = as.data.frame(scale(X, center = min, scale = max - min))

XY_n <- cbind(Y, Xn)
return(XY_n)
}

#Full libraries without scores
tetra_160n <- normalize(tetra_160)
write.csv(tetra_160n, "./data/Tetra_prop_n.csv", row.names = F)

```

## Function: PCAprop
The function PCAprop uses the sequence properties from the library _Peptides_. It takes the following arguments:

* **_n_**: Number of residues per sequence. n = 4 for tetra-peptide library, n = 5 for penta-peptide library, and so on. 
* **_pH_**: pH used to calculate the charge. DEFAULT = 7.4.
* **_custom.list_**: Boolean element indicating whether you want to use a customize library instead of the one generated with _n_. If TRUE, you need to indicate the library in PeList. DEFAULT = FALSE. 
* **_PeList_**: A vector or data.frame containing the customize library you want to obtain the properties of. Only necesary if you do _custom.list_ = *TRUE*. 
* **_rem.cys_, _rem.met_, _rem.sali_**: Boolean indicating whether you want to remove the sequences with Cysteine (_rem.cys_), Methionine (_rem.met_), or with 2 or more small aliphatic residues (_rem.sali_), or not. DEFAULT = FALSE.
* **_norm_**: Boolean indicating whether you want to normalize the properties based on min-max normalizaiton for each property. DEFAULT = F

```{r,  eval=F, PCAprop}
PCAprop <- function(n, pH = 7.4, custom.list = FALSE, PeList, rem.cys = F, rem.met = F, rem.sali = F, norm = F){
  
  system.time({
    if (custom.list == FALSE){PeList <- expand.grid(rep(list(Peptides::aaList()), n))}
    
    if (data.class(PeList) == "data.frame" & !is.null(ncol(PeList))){
      PeList <- paste(PeList[,1],
                      PeList[,2],
                      PeList[,3],
                      if (n >= 4){PeList[,4]},
                      if (n >= 5){PeList[,5]},
                      if (n >= 6){PeList[,6]},
                      if (n >= 7){PeList[,7]},
                      sep = "")
    }
    
    if (rem.cys == TRUE){PeList <- PeList[!grepl("C", PeList)]} # Remove sequences with Cys
    
    if (rem.met == TRUE){PeList <- PeList[!grepl("M", PeList)]} # Remove sequences with Met
    
    if (rem.sali == TRUE) {
      
      remove <- which(nchar(gsub("[AGILV]", "", PeList)) <= 2)  # Remove sequences with 2 or more small aliphatic amino acids
     
       if (length(remove) != 0) {
         PeList <- PeList[-remove]
       }
      }
    
    peptides <- cbind(PeList,
                      as.data.frame(do.call(rbind, Peptides::crucianiProperties(PeList))),
                      as.data.frame(do.call(rbind, Peptides::fasgaiVectors(PeList))),
                      as.data.frame(do.call(rbind, Peptides::kideraFactors(PeList))),
                      as.data.frame(do.call(rbind, Peptides::protFP(PeList))),
                      as.data.frame(do.call(rbind, Peptides::tScales(PeList))),
                      as.data.frame(do.call(rbind, Peptides::vhseScales(PeList))),
                      as.data.frame(do.call(rbind, Peptides::zScales(PeList))))
    
    pI <- cbind(pI_EMBOS = Peptides::pI(PeList, pKscale = "EMBOSS"),
                pI_Bjellqvist = Peptides::pI(PeList, pKscale = "Bjellqvist"),
                pI_Lehninger = Peptides::pI(PeList, pKscale = "Lehninger"),
                pI_Murray = Peptides::pI(PeList, pKscale = "Murray"),
                pI_Rodwell = Peptides::pI(PeList, pKscale = "Rodwell"),
                pI_Sillero = Peptides::pI(PeList, pKscale = "Sillero"))
    
    mShft_15n <- Peptides::massShift(PeList, label = "15n", aaShift = NULL, monoisotopic = TRUE)
    
    charges <- cbind(Peptides::charge(PeList, pH = pH, pKscale = "EMBOSS"))
    
    hydrophobicity_indexes <- c("Aboderin", "AbrahamLeo", "BlackMould", 
                                "BullBreese", "Casari", "Chothia", "Cid", 
                                "Cowan3.4", "Cowan7.5", "Eisenberg", "Fasman", 
                                "Fauchere", "Goldsack", "Guy", "HoppWoods", 
                                "interfaceScale_pH2", "interfaceScale_pH8", "Janin", 
                                "Jones", "Juretic", "Kuhn", "KyteDoolittle", 
                                "Levitt", "Manavalan", "Miyazawa", "octanolScale_pH2", 
                                "oiScale_pH2", "oiScale_pH8", "Parker", "Ponnuswamy", 
                                "Prabhakaran", "Rao", "Rose", "Roseman", "Sweet", 
                                "Tanford", "Welling", "Wilson", "Wolfenden", 
                                "Zimmerman")
    
    hydrophobicity <- list()
    for (i in 1:length(hydrophobicity_indexes)) {
      hydrophobicity[[i]] <- Peptides::hydrophobicity (PeList, scale = hydrophobicity_indexes[i])
      names(hydrophobicity)[[i]] <- paste0("hb_",hydrophobicity_indexes[i])
    }
    
    
    peptides <- cbind(peptides, data.frame(aIndex = Peptides::aIndex(PeList),
                                           Boman = Peptides::boman(PeList),
                                           chrg_EMBOSS = charges,
                                           hmoment1 = Peptides::hmoment(PeList, angle = 100, window = 11),
                                           hmoment2 = Peptides::hmoment(PeList, angle = 160, window = 11),
                                           hydrophobicity,
                                           instaIndex = Peptides::instaIndex(PeList),
                                           mShft_15n,
                                           mw1 = Peptides::mw(PeList, monoisotopic = TRUE),
                                           pI))
    if (norm == T) {
        X <- peptides[,-1]
        Sequence <- peptides[,1]
          
        max = apply(X , 2 , max)
        min = apply(X, 2 , min)
        Xn = as.data.frame(scale(X, center = min, scale = max - min))
        
        peptides <- cbind(Sequence, Xn)
      }
    
    return(peptides)
  })
}

```

