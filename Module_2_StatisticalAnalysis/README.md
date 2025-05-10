This work aims to use R language for systematic analysis of data. The main objectives include:
- Missing values Identification and Imputation (lipids5 as an example);
- A probability distribution was fitted to each variable by `outcome` grouping;
- Extract the parameters of the optimal distribution of each variable (e.g., mean, standard deviation, λ, etc.);
- Perform significance tests between groups (`Brunner-Munzel`);
- The consolidated results are output as a structured summary table `./Outputs/all_results.csv`;
- Use EDA reports to automatically summarize data quality and structure (`./Outputs/all_results.html`)



Highlights of processing `lipids5` variable:
- The initial missing rate is about 24%;
- The PMM method was successfully used for interpolation;
- The distribution curves almost coincide before and after interpolation → support the use of the PMM method;
- After imputation, the distribution was fitted, which was a lognormal distribution




