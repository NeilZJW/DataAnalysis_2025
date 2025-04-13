# Missing Data Imputation and Outlier Detection

This repository presents a complete workflow for handling missing values and identifying potential error records in a clinical dataset using R. The techniques include MCAR testing, multiple imputation with MICE, and outlier detection via Local Outlier Factor (LOF) and z-score diagnostics.

The working language will be R, based on the scripts provided in the course materials.

Dataset Source: Course materials: Module 1

Source code: ./DataClearning/

Output: ./Outputs/graphs/

# Working flows
1. **Data Import & Cleaning**
   - Removed irrelevant variables
   - Filtered variables with <=35% missing values
2. **Missing Data Analysis**
   - Visualized missing patterns with `visdat`, `naniar`
   - Applied `mcar_test()` to evaluate missingness mechanism
   - Concluded data is **not MCAR**
3. **Multiple Imputation**
   - Used `mice` package with:
     - `method = "rf"` (Random Forest)
     - `method = "pmm"` (Predictive Mean Matching)
   - Generated and compared imputed datasets
4.  **Visualization**
   - Plotted distributions: original vs. imputed
   -  Compared distributions of imputed vs. original data, both PMM and RF methods produced similar results (Compared each feature)
   - In this project, `PPM` is selected as the recommended interpolation method, mainly considering its ability to maintain the original data distribution. The actual effect is close to the `RF` method.
5.  **Outlier Detection**
  - Applied LOF (`dbscan::lof`) to detect multivariate anomalies
  - Visualized top LOF scores via histogram and scatter plot
  - Integrated z-score method to pinpoint specific suspicious fields
6.  **Suspected Error Reporting**
  - Extracted high-LOF records
  - Flagged extreme values per variable (z-score > 3)
  - Exported result for manual review - Future work, through screening methods, has found suspicious outliers that require human intervention to judge whether they are wrong values


# Insights
- **pmm** preserves original variable distribution better than **rf**
- LOF combined with z-score offers robust error detection
- Manual review is essential to differentiate between true outliers and data entry errors


