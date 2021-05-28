# 100_lakes_2019

This are the script and data used in "Factors of Biodegradability of Dissolved Organic Matter in Lake Water". 

Run "Scripts" in this order: 
1.respiration_rate.R : computes the respiration rate in all samples. 
  - loads the biodegradability raw data from Data/All + "All_samples.csv" (contains the ID of the samples)
  - uses "Presens_functions.R" (to clean raw data) and "doxtdt.R" (computes the respiration rate) to analyse them
  - returns "all_RR.csv" (all respiration rates with replicates), "id_table.csv" (matching IDs), "RRtotal.csv" (summary of respiration rate for a 1000 lakes) and "RR100s.csv" (summary of the respiration rate of the 73 lakes used in the article)
  
 2.prepare_data.R : merge the respiration rate data with water chemistry/DNOM data
  - loads "RR100.csv" and "CBA100Lakes_Master.xlsx"
  - returns "rr100lakes.csv"
  
 3.data_analysis.rmd : computes models and generates figures
  - loads "rr100lakes.csv"
  - returns "lakes73u.xlsx" with the selected covariates
  - returns all the figures of the article and supplementary material in .png
  - returns "3.data_analysis.html", an HTML version of the markdown
  
