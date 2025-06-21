# CHD Risk Prediction based on Framingham Dataset
This is a complete data analysis project based on the Framingham Heart Study dataset, aiming to predict the possibility of developing coronary heart disease (CHD) in the next 10 years. The project covers the entire process from data cleaning, variable selection, correlation analysis, to logistic regression modeling and model evaluation.

This project uses the publicly available dataset from the Framingham Heart Study, which contains health data of over 4,000 respondents. 
Target variable: 
- TenYearCHD: Will coronary heart disease be developed in the next ten years? (1 = Yes, 0 = No) 
- Input variables include multiple indicators such as age, gender, blood pressure, diabetes, BMI, cholesterol, smoking, and blood sugar levels.

# Analysis process
## 1. Data processing and imputation
- Use the `mice` package to perform `PMM` imputation for missing values. 
- Remove or correct logical outliers (such as the daily cigarette count of non-smokers). 
## 2. Correlation Analysis
- Spearman correlation + permutation test was used to evaluate the strength of the association between independent variables and `TenYearCHD`. 
## 3. Modeling and Comparison
- Establish two models: 
  - Full-variable logistic regression model 
  - Stepwise regression model (based on AIC) 
- Variable significance test and OR (Odds Ratio) plotting analysis. 
## 4. Model Evaluation
- Compare model performance using ROC curves 
- Output the AUC value and confusion matrix. 
- Compare the model accuracy, complexity and the differences in TP and FP.

# Conclusion
## 1. Variables significantly associated with CHD risk
Through Spearman correlation test and Logistic regression analysis, we found that the following variables have statistically significant relationships with the risk of CHD: 
- Age: The risk of CHD increases significantly with age (OR = 1.07, p < 0.001).
- Systolic blood pressure (sysBP): For every 1-unit increase in systolic blood pressure, the risk of CHD increases by approximately 1% (OR = 1.01, p < 0.001).
- A history of prevalent hypertension (PrevalentHyp): Individuals with a history of hypertension have a higher risk of CHD (OR = 1.27).
- The risk of CHD in men is 1.66 times that in women (OR = 1.66, p < 0.001).
- A history of stroke (prevalent stroke): Individuals with a history of stroke have a significantly increased risk of CHD (OR = 2.76, p < 0.05).
- Glucose and daily cigarette consumption (cigsPerDay) were also positively correlated with CHD (p < 0.001).
## 2. Analysis of the Effect of LogisticRegression Model
![1750526470236](https://github.com/user-attachments/assets/e8e3c0f3-9c64-4523-8578-4511b6461016)
- Although the AIC of the stepwise regression model is lower (indicating better model quality), its ability to identify CHD patients is slightly lower (with a decrease of 5 true positives). 
- The AUC difference between the two is extremely small (<0.001), indicating comparable performance.
## 3. Visualization and Interpretation of Models
- We plotted the ROC curve and verified that the model has a good discrimination ability (AUC ≈ 0.73). 
- The Odds Ratio graph was used to help interpret the intensity and direction of the influence of each variable on CHD, enhancing the interpretability of the results. 
- The confusion matrix reveals the distribution of TP, FP, FN, and TN in the model's classification in reality, facilitating subsequent model optimization or deployment.
