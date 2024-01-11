A project titled "Integrating multiple lines of evidence to assess the effects of maternal BMI on pregnancy and perinatal outcomes". 

This repository includes the analysis plan and all code related to the manuscript.

This repository contains:
1) Directory: 1_main_analyses
   
Scripts

•	0_pool_MV.R
Pools MV regression results from each study in both a fixed and random effects meta-analysis

•	1_pool_snp-trait.R
Selects SNPs as instruments for maternal BMI and extracts corresponding exposure GWAS summary data

•	2_mr.R

2)	Directory: 2_output_plots_tables
   
Scripts:

•	0_mr_plots.R
Creates plots for MR analyses

•	1_main_plots.R
Creates the main plots showing MR x MVR and MVR x PNC analyses

•	2_supp_plots_info.R
Creates supplementary plots for MR x MVR and MVR x PNC analyses with the continuous outcomes

•	3_supp_plots_facets.R
Creates supplementary plots including by study, comparison of unadjusted/adjusted, fixed/random effects

•	4a_sens_PNC_ALSPAC.R
Creates a study specific plot comparing maternal and paternal BMI estimates for different models in the PNC analysis (similarly for scripts 4b and 4c)

•	4b_sens_PNC_GenR.R

•	4c_sens_PNC_MoBa.R

•	5_supp_tables.R
Creates a supplementary table with between study heterogeneity information by outcome

•	6_pool_sens_MVR.R
Pools data for study by pre/during pregnancy weight

•	7_plot_sens_MVR.R
Plots results by pre/during pregnancy weight
