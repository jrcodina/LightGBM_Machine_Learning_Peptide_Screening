# LightGBM Machine Learning Framework for Screening Large Peptide Libraries
Scripts for the paper *Accelerating the Screening of Small Peptide Ligands by Combining Peptide-Protein Docking and Machine Learning*. Available in: 
https://www.mdpi.com/1422-0067/24/15/12144

The data sets can either be downloaded in here and then decompressed, or by going to the publication's supplementary material section under the "Scripts_reproducibility" folder.


## Cite
Please cite this paper if you found it useful.Thank you. 


### AMA Style

Codina J-R, Mascini M, Dikici E, Deo SK, Daunert S. Accelerating the Screening of Small Peptide Ligands by Combining Peptide-Protein Docking and Machine Learning. International Journal of Molecular Sciences. 2023; 24(15):12144. https://doi.org/10.3390/ijms241512144


### BibTeX

@Article{ijms241512144,
AUTHOR = {Codina, Josep-Ramon and Mascini, Marcello and Dikici, Emre and Deo, Sapna K. and Daunert, Sylvia},
TITLE = {Accelerating the Screening of Small Peptide Ligands by Combining Peptide-Protein Docking and Machine Learning},
JOURNAL = {International Journal of Molecular Sciences},
VOLUME = {24},
YEAR = {2023},
NUMBER = {15},
ARTICLE-NUMBER = {12144},
URL = {https://www.mdpi.com/1422-0067/24/15/12144},
ISSN = {1422-0067},
ABSTRACT = {This research introduces a novel pipeline that couples machine learning (ML), and molecular docking for accelerating the process of small peptide ligand screening through the prediction of peptide-protein docking. Eight ML algorithms were analyzed for their potential. Notably, Light Gradient Boosting Machine (LightGBM), despite having comparable F1-score and accuracy to its counterparts, showcased superior computational efficiency. LightGBM was used to classify peptide-protein docking performance of the entire tetrapeptide library of 160,000 peptide ligands against four viral envelope proteins. The library was classified into two groups, &lsquo;better performers&rsquo; and &lsquo;worse performers&rsquo;. By training the LightGBM algorithm on just 1% of the tetrapeptide library, we successfully classified the remaining 99%with an accuracy range of 0.81&ndash;0.85 and an F1-score between 0.58&ndash;0.67. Three different molecular docking software were used to prove that the process is not software dependent. With an adjustable probability threshold (from 0.5 to 0.95), the process could be accelerated by a factor of at least 10-fold and still get 90&ndash;95% concurrence with the method without ML. This study validates the efficiency of machine learning coupled to molecular docking in rapidly identifying top peptides without relying on high-performance computing power, making it an effective tool for screening potential bioactive compounds.},
DOI = {10.3390/ijms241512144}
}

