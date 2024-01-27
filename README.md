# IPCA-thesis
These codes are for Tony Wu's Mannaheim MMM master's thesis. 

1. characterisrtics_data_Sep2023.m

   * generates data from 'characteristics_data_feb2017.m' downloaded from Seth Pruitt's Dropbox: [https://www.dropbox.com/sh/sykn8aac1e2s6ry/AAAsGmD_qKRqDpOQhRAyWir_a?dl=0]
   
2. datamaker.m, IPCA_empirical_datacall_FNW36.m

   * datamaker.m prepares data for further analysis, it can be customized in the option panel. IPCA_empirical_datacall_FNW36 is served as a functon in datamaker

3. IPCA_empirical_GB.m,IPCA_empirical_GBGA.m ,stocklevelpca.m ,managedportPCA.m

   * IPCA_empirical_GB.m,IPCA_empirical_GBGA.m are use for analyzing the performance of frestricted and unrestricted fundamental IPCA model, respectively.
   * stocklevelpca.m ,managedportPCA.m are for analysis of PCA model at stock and managed-portfolio level

4. Observable_Factor_Regression.m, IPCA_empirical_GBGDFF.m, IPCA_empirical_OtherFactors

   * Observable_Factor_Regression.m processes the analysis of FF3, FFC4, FF5, FF6, BY, HXZ, M4 models  without instrument
   * IPCA_empirical_GBGDFF.m processes the analysis of restricted fundamental model,restricted fundamental model+FFC 6 factors, and instrumented FFC 6 factors
   * IPCA_empirical_OtherFactors.m processes the analysis of restricted fundamental model,restricted fundamental model+other obserable factors, and instrumented observable factors  

5. fred-database_code/ fredfactorprocess.m, Welchcode.m, LASSOselection.m

   * fredfactorprocess.m, Welchcode.m prepare the stationary macroeconomics features, including FRED-MD data '2015-01.csv / 2023-09.csv' , Welch data '2020Welch', and median of firm characteristics
   * LASSOselection.m applies LASSO regression for macroeconomic factor selection, the selected features are saved under 'macro_1.csv'

6. fred-database_code/ RNNLSTM.ipynb, modelfit.ipynb

   * RNNLSTM.ipynb constructs Recurrent Neural Network with LSTM cell. The time steps could be modified under 'stepin, stepout'. It provides the extracted hidden state from last timestep and plots the comparison between the extracted hidden state to NBER recession date
   * modelfit.ipynb plots the contribution of each macroeconoic features within the neural network

7. IPCA_empirical_puremacro.m, IPCA_empiricla_GBGDmacro.m

   * IPCA_empirical_puremacro. reports the performance of macroeconomic features selected by LASSO regression w/o passing RNN
   * IPCA_empiricla_GBGDmacro.m reports the performance of the last hidden state

8. IPCA_empirical_GBGDFFtest.m, IPCA_empirical_GBGDmacrotest.m
   
   *  IPCA_empirical_GBGDFFtest.m, IPCA_empirical_GBGDmacrotest.m conduct the Wald-type test for FF factors, macroeconomic factors from last hidden state, macroeconomic factors from LASSO selection, respectively.

9. Plotloadings.m, Plotmacroloadings.m
    
   * Plotloadings.m and Plotmacroloadings.m plot the weight of firm characteristics 
    for fundamental IPCA and extended IPCA, respectively

10. IPCA_empirical_GB_outofsample.m, ,IPCA_empirical_GBGDmacro_outofsample.m, tanpf.m, tanpfnext.m

  * IPCA_empirical_GB_outofsample.m, ,IPCA_empirical_GBGDmacro_outofsample.m generated the out-of-sample performance (fits and Sharpe ratio) of fundamental and extended model
  * tanpf.m, tanpfnext.m are served as function for outofsample estimation

11. Table_OutofSample_Fits.m, Table_OutofSamplemacro_Fits.m

  * Table_OutofSample_Fits.m, Table_OutofSamplemacro_Fits.m generate tables for out-of-sample performance of fundamental, extended IPCA, and FFC6 model

12. individualtest.m, individualtest_macro.m 

  *  individualtest.m, individualtest_macro.m serve as the Wald-type test for individual instruments for fundamental IPCA and extended IPCA, respectively





