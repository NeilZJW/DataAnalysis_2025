# Statistical Analysis 2

This work is a complete statistical analysis of the `lipids5` variable using R, completed as part of a practical class assignment in **Statistical Analysis**.

## Analysis Workflow

### 1. Normality Testing
- Method: **Shapiro-Wilk Test**
- Conclusion: Most variables are **non-normally distributed** (p < 0.05)

### 2. Correlation Analysis
- Method: **Spearman correlation** with **Permutation Test** for significance
- Output: Correlation coefficient table + p-values
- Significant variables (p < 0.05):  
  `lipids1`, `lipids2`, `lipids3`, `carb_metabolism`, `hormone2`, `lipid_pero4`
- Note: `carb_metabolism`, c is not an English letter, that is, a Russian с

### 3. Regression Analysis
- Models Tested:
  - Linear Regression
  - Polynomial (2nd and 3rd degree)
  - Exponential
  - Logarithmic
- Model selection criterion: **BIC (Bayesian Information Criterion)**
- Output: Best model per variable based on lowest BIC

### 4. Visualization
- Tool: Base R plotting system (`plot()`, `lines()`, `abline()`)
- Each significant variable plotted against `lipids5` with:
  - Black points (raw data)
  - Blue lines (data trends)
  - Red regression line
