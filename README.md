# IPCA-thesis
These codes are for Tony Wu's Mannaheim MMM master's thesis. 

1. characterisrtics_data_Sep2023.m

   * generates data from 'characteristics_data_feb2017.m' downloaded from Seth Pruitt's Dropbox: [https://www.dropbox.com/sh/sykn8aac1e2s6ry/AAAsGmD_qKRqDpOQhRAyWir_a?dl=0]
   
2. datamaker.m, IPCA_empirical_datacall_FNW36.m

   * datamaker.m prepares data for further analysis, it can be customized in the option panel. IPCA_empirical_datacall_FNW36 is served as a functon in datamaker

3. IPCA_empirical_GB.m,IPCA_empirical_GBGA.m ,stocklevelpca.m ,managedportPCA.m

   * IPCA_empirical_GB.m,IPCA_empirical_GBGA.m are use for analyzing the performance of frestricted and unrestricted fundamental IPCA model, respectively.
   * stocklevelpca.m ,managedportPCA.m are for analysis of PCA model at stock and managed-portfolio level

4. Observable_Factor_Regression.m, IPCA_empirical_GBGDFF.m

   * Observable_Factor_Regression.m processes the analysis of FF3, FFC4, FF5, FF6, BY, HXZ, M4 models  without instrument
   * IPCA_empirical_GBGDFF.m processes the analysis of restricted fundamental model,restricted fundamental model+FFC 6 factors, and instrumented FFC 6 factors  

5. fred-database_code/ fredfactorprocess.m, Welchcode.m, LASSOselection.m

   * fredfactorprocess.m, Welchcode.m prepare the stationary macroeconomics features, including FRED-MD data '2015-01.csv / 2023-09.csv' , Welch data '2020Welch', and median of firm characteristics
   * LASSOselection.m applies LASSO regression for macroeconomic factor selection, the selected features are saved under 'macro_1.csv'

6. fred-database_code/ RNNLSTM.ipynb, modelfit.ipynb

   * RNNLSTM.ipynb constructs Recurrent Neural Network with LSTM cell. The time steps could be modified under 'stepin, stepout'. It provides the extracted hidden state from last timestep and plots the comparison between the extracted hidden state to NBER recession date
   * modelfit.ipynb plots the contribution of each macroeconoic features within the neural network

