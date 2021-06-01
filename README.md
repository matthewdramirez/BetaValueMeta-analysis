# BetaValueMeta-analysis

This repository contains beta (β) value datasets and R code associated with the article "Meta-analysis of primary producer amino acid δ<sup>15</sup>N values and their influence on trophic position estimation," currently in review at Methods in Ecology and Evolution. 

_Reference:_ Ramirez MD, Besser AC, Newsome SD, McMahon KW (in review) Meta-analysis of primary producer amino acid δ<sup>15</sup>N values and their influence on trophic position estimation. Methods in Ecology and Evolution XX, XXX-XXX. https://doi.org/

### File Descriptions ###

* ***Beta_data.xlsx*** includes β value data, the difference between trophic and source amino acid (AA) δ<sup>15</sup>N values within primary producers, resulting from a meta-analysis of the published primary producer literature. See _Metadata_ tab within this file for full description of the dataset. 

  Note: As of June 1, 2021, "Beta_data.xlsx" excludes unpublished data from Chen et al. 2020 (ice algae, n = 6) and A. C. Besser [arid habitat primary producer; freshwater eukaryotic microalgae (n = 4), grass (n = 7), forb (n = 6), cactus (n = 5), shrub (n = 6)]. Please direct data inquiries to Shaomin Chen (ice algae; shaomin.chen[at]dal.ca) or Alexi C. Besser (arid primary producer; acbesser[at]unm.edu). 

* ***Fig1 - Conceptual Simulation.R*** includes code to reproduce Figure 1B, which is a simulation demonstrating how variation in β values and AA-specific TDFs propagate through a simple hypothetical food chain to influence consumer trophic position estimates.

* ***Fig5 - Sensitivity Analysis.R*** includes code to reproduce Figure 5, which illustrates how consumer trophic position estimates change as a function of variation in mean β values and trophic discrimination factors. 

  Note: This code produces consumer-specific results panels that were compiled by the authors using Adobe Illustrator.


