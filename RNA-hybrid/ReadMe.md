# RNA-hybrid Circuit Folders:

(1) Sensitivity Analysis

This folder containes only two files: Protein_Hill_Model.m and Protein_Sensitivity_Analysis.m. The Protein_Hill_Model.m file contains the ODEs used to model the RNA-hybrid Circuit. The Sensitivity_Analysis.m begins by defininng the nominal values for the respective parameters. The first section after the varible definition "Nominal Value Calculation" generates the nominal output values. 

The following section titled "Sensativity Analysis" sets a range +/- 50% of the nominal value from which a value will be randomly selected using the function randi. The range of inducer inputs uitlized for the analysis are as follows: aTc  (200, 100, 20, 2) ng/ml and IPTG (1, 0.4, 0.02, 0.01) nM. 
		
(2) Protein Fitting		

This folder containes files used for two objectives: 1) Train the model using experimental data 2) Visualization/Experimental Validation. The training or fitting of the model was done using three files: Control_Model_Fitting.m, Obj.m, and Protein_Hill_Model.m. As mentioned previously, the Protein_Hill_Model.m contains the ODEs used to describe the circuit. The Obj.m file contains the objective function, which calculates the difference or error between the simulated output of the model and the experimental values. These experimental values are stored in a file titled Exp_Data_V4.mat,  which contains a 4x4 training set. The total run time of the experiments was 480 minutes; however, the first 50 min were excluded in the data fitting. 
		
The Control_Model_Fitting.m file contains the inital guess values for each respective parameter. Additionaly, "fmincon" was uitlized to train the model by calling the Obj.m file. Below the "fmincon" function are the defined upper and lower bounds for the parameters. 

The results from the fitting were saved as a .mat file titled Final_Fitting_Results.mat and then used in the Protein_Fitting_Visulization.m file to plot and compare the model simulation output with the training and additonal experimental data which is stored in Exp_Valid_V1.mat.  
		
		
		
		
		
