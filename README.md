# RNA-Only Circuit Model

The are two files for simulation: The main file named RNA_ONLY_main.m and the ODE file named RNA_ONLY_ODE.m
All the pre-processed experimental data (GFP measurement) are stored in the provided .mat files.
For RNA-only circuit, the training data are the insulated Y data, which have been collated into the variable ExpData
Each row in ExpData represents the GFP measurement for different IPTG and Arabinose concentrations.
The lookup value below shows the associated IPTG and Arabinose concentrations for each row in ExpData
Note that a '1' for Arabinose represents 6600mM and a '1' for IPTG represents 1000mM

Row 1: Arabinose = 1, IPTG = 1/2, Row 2: Arabinose = 1, IPTG = 1/8, Row 3: Arabinose = 1, IPTG = 1/32, Row 4: Arabinose = 1, IPTG = 1/128
Row 5: Arabinose = 1/4, IPTG = 1/2, Row 6: Arabinose = 1/4, IPTG = 1/8, Row 7: Arabinose = 1/4, IPTG = 1/32, Row 8: Arabinose = 1/4, IPTG = 1/128
Row 9: Arabinose = 1/16, IPTG = 1/2, Row 10: Arabinose = 1/16, IPTG = 1/8, Row 11: Arabinose = 1/16, IPTG = 1/32, Row 12: Arabinose = 1/16, IPTG = 1/128
Row 13: Arabinose = 1/64, IPTG = 1/2, Row 14: Arabinose = 1/64, IPTG = 1/8, Row 15: Arabinose = 1/64, IPTG = 1/32, Row 16: Arabinose = 1/64, IPTG = 1/128

Also, as the first 90 minutes data are not considered in the training due to low OD600 value, we data used for training begins from row 7. 
For example, to use data set with Arabinose = 1/4 (1650mM) and IPTG = 1/8 (125mM), we denote this as ExpData(6,7:end)
