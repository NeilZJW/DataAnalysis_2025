#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset--------------------------
example_df<-read.csv("distribution.csv", header = TRUE,dec = ',', sep = ";")
factor_df <- read.csv("factor_data.csv")
imputed_df <- read.csv("imputed_data.csv")
# Display structure with variable types
str(example_df)
str(factor_df)
str(imputed_df)
#---------------merge two files-------------------------
data_for_analysis <- merge(
  factor_df, 
  imputed_df, 
  by = "record_id",        # column for merge
  all = FALSE       # FALSE = INNER JOIN (only coincidences), TRUE = FULL JOIN
)
str(data_for_analysis)

# save data_for_analysis in CSV
write.csv(data_for_analysis, "data_for_analysis.csv", row.names = FALSE)  

#------------------Probability Distributions----------------------- 
install.packages("MASS", dependencies=T)
library(MASS)
#----------------example-----------------------------------------
summary (example_df)
example_df$value <- as.numeric(example_df$value)
summary(example_df)
#building histograms for example
# normal distribution
val<-example_df[example_df$distribution=="norm",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="normal")

fit
#lognormal distribution
val<-example_df[example_df$distribution=="lognorm",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="lognormal")

fit

unname(fit$estimate[1])

unname(fit$estimate[2])

m_log<-exp(unname(fit$estimate[1]))*sqrt(exp(unname(fit$estimate[2])^2))
m_log
sd_log<-sqrt(exp(2*unname(fit$estimate[1]))*(exp(unname(fit$estimate[2])^2)-1)*sqrt(exp(unname(fit$estimate[2])^2)))
sd_log
#exponential distribution

val<-example_df[example_df$distribution=="exp",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="exponential")

fit

unname(fit$estimate[1])

m_exp<-1/unname(fit$estimate[1])
m_exp
#Poisson distribution
val<-example_df[example_df$distribution=="pois",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="Poisson")

fit

unname(fit$estimate[1])

sd_pois<-sqrt(unname(fit$estimate[1]))
sd_pois
#Selecting a Distribution Model

val<-example_df[example_df$distribution=="lognorm",]$value

fit_1<-fitdistr(val, densfun="normal")
fit_2<-fitdistr(val, densfun="lognormal")
fit_3<-fitdistr(val, densfun="exponential")

#Bayesian Information Criterion calculation
BIC(fit_3)

#calculation of the Bayesian information criterion for all models
BIC_value<-c(BIC(fit_1), BIC(fit_2), BIC(fit_3))

#forming a vector with the name of the models
distribution<-c("normal", "lognormal", "exponential")

#combining the results into a final table
rez<-data.frame(BIC_value=BIC_value, distribution=distribution)

#sort table in ascending order of Bayesian Information Criterion value
rez<-rez[order(rez$BIC_value, decreasing=F),]

rez


#calculation of absolute values of the confidence interval for the mean of a lognormal distribution
error_min<-unname(fit_2$estimate[1])-unname(fit_2$sd[1])
error_max<-unname(fit_2$estimate[1])+unname(fit_2$sd[1])

error_min
error_max

m<-exp(unname(fit_2$estimate[1]))*sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_min<-exp(error_min)*sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_max<-exp(error_max)*sqrt(exp(unname(fit_2$estimate[2])^2))

value_error_min
m
value_error_max

#--------------data for analysis--------------------------
#building histograms
value_d1<-data_for_analysis$lipids1
hist(value_d1)
value_d2<-data_for_analysis$lipids2
hist(value_d2)
value_d3<-data_for_analysis$lipids3
hist(value_d3)
value_d4<-data_for_analysis$lipids4
hist(value_d4)


# d1 distribution estimate


fit_d1_1<-fitdistr(value_d1,densfun="normal")
fit_d1_2<-fitdistr(value_d1,densfun="lognormal")
fit_d1_3<-fitdistr(value_d1,densfun="exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d1

BIC_value_d1 <- c(BIC(fit_d1_1),BIC(fit_d1_2),BIC(fit_d1_3))
distribution <-c("normal","lognormal","exponential")
result_d1<-data.frame(BIC_value_d1=BIC_value_d1, distribution=distribution)
result_d1
min(result_d1$BIC_value_d1)
distribution_d1<-result_d1[result_d1$BIC_value_d1==min(result_d1$BIC_value_d1),]$distribution
distribution_d1
# Finding parameters for d1
fit_d1_1$estimate[1:2]

# d2 distribution estimate


fit_d2_1<-fitdistr(value_d2,densfun="normal")
fit_d2_2<-fitdistr(value_d2,densfun="lognormal")
fit_d2_3<-fitdistr(value_d2,densfun="exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d2

BIC_value_d2 <- c(BIC(fit_d2_1),BIC(fit_d2_2),BIC(fit_d2_3))
distribution <-c("normal","lognormal","exponential")
result_d2<-data.frame(BIC_value_d2=BIC_value_d2, distribution=distribution)
result_d2
min(result_d2$BIC_value_d2)
distribution_d2<-result_d2[result_d2$BIC_value_d2==min(result_d2$BIC_value_d2),]$distribution
distribution_d2
# Finding parameters for d2
fit_d2_1$estimate[1:2]

#-----------descriptive statistics------------------
#-----------for publication tables-----------------
install.packages("gtsummary")
library(gtsummary)

tbl_summary(data_for_analysis)  # Automatic table
tbl_summary(data_for_analysis, by = outcome)  # By groups
#---------------Creating a custom table--------------
# Homework: Creating a custom table with descriptive statistics results
# 加载必要包
library(MASS)
library(lawstat)

# 封装：分析单个变量函数
analyze_variable <- function(varname, data) {
    # 提取并清洗两组
    val_0 <- data[[varname]][data$outcome == "0"]
    val_1 <- data[[varname]][data$outcome == "1"]
    val_0 <- val_0[is.finite(val_0) & !is.na(val_0)]
    val_1 <- val_1[is.finite(val_1) & !is.na(val_1)]
    
    # 样本数不足时跳过
    if (length(val_0) < 10 || length(val_1) < 10) {
        return(NULL)
    }
    
    # 拟合 outcome == 0
    fits_0 <- list()
    fits_0$normal <- fitdistr(val_0, "normal")
    fits_0$exponential <- fitdistr(val_0, "exponential")
    if (all(val_0 > 0)) {
        fits_0$lognormal <- fitdistr(val_0, "lognormal")
    }
    
    # BIC比较 outcome == 0
    bics_0 <- sapply(fits_0, BIC)
    dist_0 <- names(which.min(bics_0))
    est_0 <- unname(fits_0[[dist_0]]$estimate)
    
    # 拟合 outcome == 1
    fits_1 <- list()
    fits_1$normal <- fitdistr(val_1, "normal")
    fits_1$exponential <- fitdistr(val_1, "exponential")
    if (all(val_1 > 0)) {
        fits_1$lognormal <- fitdistr(val_1, "lognormal")
    }
    
    # BIC比较 outcome == 1
    bics_1 <- sapply(fits_1, BIC)
    dist_1 <- names(which.min(bics_1))
    est_1 <- unname(fits_1[[dist_1]]$estimate)
    
    # Brunner-Munzel检验
    pval <- brunner.munzel.test(val_0, val_1)$p.value
    
    # 确定 outcome = 0 的参数名
    if (dist_0 == "normal") {
        param1_name_0 <- "Mean_0"
        param2_name_0 <- "SD_0"
    } else if (dist_0 == "lognormal") {
        param1_name_0 <- "Mu_0"
        param2_name_0 <- "Sigma_0"
    } else if (dist_0 == "exponential") {
        param1_name_0 <- "Rate_0"
        param2_name_0 <- NA
    }
    
    # 确定 outcome = 1 的参数名
    if (dist_1 == "normal") {
        param1_name_1 <- "Mean_1"
        param2_name_1 <- "SD_1"
    } else if (dist_1 == "lognormal") {
        param1_name_1 <- "Mu_1"
        param2_name_1 <- "Sigma_1"
    } else if (dist_1 == "exponential") {
        param1_name_1 <- "Rate_1"
        param2_name_1 <- NA
    }
    
    # 构建结果表格
    result <- data.frame(
        Variable = varname,
        Outcome_0_Distribution = dist_0,
        Outcome_0_ParamName1 = param1_name_0,
        Outcome_0_Param1 = round(est_0[1], 4),
        Outcome_0_ParamName2 = param2_name_0,
        Outcome_0_Param2 = ifelse(length(est_0) > 1, round(est_0[2], 4), NA),
        Outcome_1_Distribution = dist_1,
        Outcome_1_ParamName1 = param1_name_1,
        Outcome_1_Param1 = round(est_1[1], 4),
        Outcome_1_ParamName2 = param2_name_1,
        Outcome_1_Param2 = ifelse(length(est_1) > 1, round(est_1[2], 4), NA),
        Brunner_Munzel_P = round(pval, 5)
    )
    return(result)
}
variables <- c("lipids1", "lipids2", "lipids3", "lipids4", 
               "hormone1", "hormone2", "hormone3", "hormone4")
all_results <- do.call(rbind, lapply(variables, analyze_variable, data = data_for_analysis))
View(all_results)  # 或 print(all_results)
library(gtsummary)
tbl_summary(all_results, by="Variable")

library(DataExplorer)
create_report(all_results)  # Generates HTML report with graphs and statistics
create_report(
    data = all_results,
    output_file = "all_results.html",  
    output_dir = getwd(),                
    report_title = "all_results"          
)

#--------------Statistical Tests---------------------
value_outcome1<-data_for_analysis[data_for_analysis$outcome=="1",]$lipids1
hist(value_outcome1)
value_outcome0<-data_for_analysis[data_for_analysis$outcome=="0",]$lipids1
hist(value_outcome0)
#-------Levene's Test for Homogeneity of Variance--------------
install.packages("car")
library(car)
str(data_for_analysis)
data_for_analysis$outcome<- as.factor(data_for_analysis$outcome)
car::leveneTest(lipids1 ~ outcome, data = data_for_analysis)
#---------------Application of the Brunner-Munzel test----------
install.packages("lawstat")
library(lawstat)
group1 <- data_for_analysis$lipids1[data_for_analysis$outcome == "0"]
group2 <- data_for_analysis$lipids1[data_for_analysis$outcome == "1"]

brunner.munzel.test(group1, group2)
#-------------comparison of results with other tests--------------
t.test(group1, group2)
wilcox.test(group1, group2)

#----------------------------EDA----------------------------------
install.packages("DataExplorer")
library(DataExplorer)
create_report(data_for_analysis)  # Generates HTML report with graphs and statistics
create_report(
  data = data_for_analysis,
  output_file = "EDA_Report.html",  
  output_dir = getwd(),                
  report_title = "EDA Report"          
)

# Homework: extra
missing_ratio <- colMeans(is.na(data_for_analysis))
missing_ratio[order(-missing_ratio)]  # 从高到低排序
# factor_pcos: 85%
data_for_analysis <- data_for_analysis[, !(names(data_for_analysis) %in% "factor_pcos")]

summary(data_for_analysis$lipids5)
sum(is.na(data_for_analysis$lipids5))
data_backup <- data_for_analysis

data_for_analysis$lipids5 <- as.numeric(data_for_analysis$lipids5)
data_for_analysis$lipids5[is.nan(data_for_analysis$lipids5)] <- NA
data_for_analysis$lipids5[is.infinite(data_for_analysis$lipids5)] <- NA



# ---- MCAR 检验 ----
install.packages("mice")
install.packages("naniar")
install.packages("dplyr")
install.packages("ggplot2")
library(mice)
library(naniar)
library(dplyr)
library(ggplot2)

# 使用的数据框
handle_MD_df <- data_for_analysis

# MCAR 检验
mcar_result <- mcar_test(handle_MD_df)

interpret_mcar <- function(mcar_result) {
    p <- mcar_result$p.value
    if (p > 0.05) {
        message("- value > 0.05 → Data is likely MCAR. Safe to delete or impute.")
    } else {
        message("- value <= 0.05 → Data is NOT MCAR. Assume MAR or MNAR.")
        message("- 建议使用多重插补（mice）方法，如 pmm / rf 等。")
    }
}
interpret_mcar(mcar_result)

# ---- 多重插补 ----
# 方法一：Random Forest
# 初始化 mice，提取默认 predictorMatrix
ini <- mice(data_for_analysis, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method

# 设置方法，只插补 lipids5
meth[] <- ""  # 默认不插补任何列
meth["lipids5"] <- "rf"

# 设置预测变量：
# 这里让 lipids5 被以下变量预测
pred["lipids5", ] <- 0
pred["lipids5", c("lipids1", "lipids2", "lipids3", "lipids4", "hormone1", "hormone2", "hormone3", "hormone4")] <- 1

# 执行插补
imp_rf <- mice(data_for_analysis, method = meth, predictorMatrix = pred, m = 1, seed = 123)
data_rf <- complete(imp_rf)

# 方法二：Predictive Mean Matching (PMM)
# 初始化 mice，提取默认 predictorMatrix
ini <- mice(data_for_analysis, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method
# 设置方法，只插补 lipids5
meth[] <- ""  # 默认不插补任何列
meth["lipids5"] <- "pmm"

# 设置预测变量：
# 这里让 lipids5 被以下变量预测
pred["lipids5", ] <- 0
pred["lipids5", c("lipids1", "lipids2", "lipids3", "lipids4", "hormone1", "hormone2", "hormone3", "hormone4")] <- 1

imp_pmm <- mice(data_for_analysis, method = "pmm", m = 1, seed = 123)
data_pmm <- complete(imp_pmm)


plot_imputation_comparison <- function(var_name) {
    original <- data_for_analysis %>%
        select(all_of(var_name)) %>%
        mutate(source = "Original")
    
    rf_imputed <- data_rf %>%
        select(all_of(var_name)) %>%
        mutate(source = "Random Forest")
    
    pmm_imputed <- data_pmm %>%
        select(all_of(var_name)) %>%
        mutate(source = "PMM")
    
    combined <- bind_rows(original, rf_imputed, pmm_imputed)
    
    ggplot(combined, aes_string(x = var_name, fill = "source", color = "source")) +
        geom_density(alpha = 0.4, size = 1, na.rm = TRUE) +
        labs(title = paste("Density Comparison of:", var_name),
             x = var_name, y = "Density") +
        theme_minimal()
}

# 运行对比图
plot_imputation_comparison("lipids5")

library(mice)

# 确保类型正确
data_for_analysis$lipids5 <- as.numeric(data_for_analysis$lipids5)
data_for_analysis$lipids5[is.nan(data_for_analysis$lipids5)] <- NA
data_for_analysis$lipids5[is.infinite(data_for_analysis$lipids5)] <- NA

# 初始化 mice，提取默认 predictorMatrix
ini <- mice(data_for_analysis, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method

# 设置方法，只插补 lipids5
meth[] <- ""  # 默认不插补任何列
meth["lipids5"] <- "pmm"

# 设置预测变量：
# 这里让 lipids5 被以下变量预测（可以再加上信任的）
pred["lipids5", ] <- 0
pred["lipids5", c("lipids1", "lipids2", "lipids3", "lipids4", "hormone1", "hormone2", "hormone3", "hormone4")] <- 1

# 执行插补
imp <- mice(data_for_analysis, method = meth, predictorMatrix = pred, m = 1, seed = 123)
data_imputed <- complete(imp)

# 替换 lipids5
data_for_analysis$lipids5 <- data_imputed$lipids5

# 验证是否成功
sum(is.na(data_for_analysis$lipids5))  # 应该是 0

# 加载必要包
library(MASS)
library(lawstat)

# 封装：分析单个变量函数
analyze_variable <- function(varname, data) {
    # 提取并清洗两组
    val_0 <- data[[varname]][data$outcome == "0"]
    val_1 <- data[[varname]][data$outcome == "1"]
    val_0 <- val_0[is.finite(val_0) & !is.na(val_0)]
    val_1 <- val_1[is.finite(val_1) & !is.na(val_1)]
    
    # 样本数不足时跳过
    if (length(val_0) < 10 || length(val_1) < 10) {
        return(NULL)
    }
    
    # 拟合 outcome == 0
    fits_0 <- list()
    fits_0$normal <- fitdistr(val_0, "normal")
    fits_0$exponential <- fitdistr(val_0, "exponential")
    if (all(val_0 > 0)) {
        fits_0$lognormal <- fitdistr(val_0, "lognormal")
    }
    
    # BIC比较 outcome == 0
    bics_0 <- sapply(fits_0, BIC)
    dist_0 <- names(which.min(bics_0))
    est_0 <- unname(fits_0[[dist_0]]$estimate)
    
    # 拟合 outcome == 1
    fits_1 <- list()
    fits_1$normal <- fitdistr(val_1, "normal")
    fits_1$exponential <- fitdistr(val_1, "exponential")
    if (all(val_1 > 0)) {
        fits_1$lognormal <- fitdistr(val_1, "lognormal")
    }
    
    # BIC比较 outcome == 1
    bics_1 <- sapply(fits_1, BIC)
    dist_1 <- names(which.min(bics_1))
    est_1 <- unname(fits_1[[dist_1]]$estimate)
    
    # Brunner-Munzel检验
    pval <- brunner.munzel.test(val_0, val_1)$p.value
    
    # 确定 outcome = 0 的参数名
    if (dist_0 == "normal") {
        param1_name_0 <- "Mean_0"
        param2_name_0 <- "SD_0"
    } else if (dist_0 == "lognormal") {
        param1_name_0 <- "Mu_0"
        param2_name_0 <- "Sigma_0"
    } else if (dist_0 == "exponential") {
        param1_name_0 <- "Rate_0"
        param2_name_0 <- NA
    }
    
    # 确定 outcome = 1 的参数名
    if (dist_1 == "normal") {
        param1_name_1 <- "Mean_1"
        param2_name_1 <- "SD_1"
    } else if (dist_1 == "lognormal") {
        param1_name_1 <- "Mu_1"
        param2_name_1 <- "Sigma_1"
    } else if (dist_1 == "exponential") {
        param1_name_1 <- "Rate_1"
        param2_name_1 <- NA
    }
    
    # 构建结果表格
    result <- data.frame(
        Variable = varname,
        Outcome_0_Distribution = dist_0,
        Outcome_0_ParamName1 = param1_name_0,
        Outcome_0_Param1 = round(est_0[1], 4),
        Outcome_0_ParamName2 = param2_name_0,
        Outcome_0_Param2 = ifelse(length(est_0) > 1, round(est_0[2], 4), NA),
        Outcome_1_Distribution = dist_1,
        Outcome_1_ParamName1 = param1_name_1,
        Outcome_1_Param1 = round(est_1[1], 4),
        Outcome_1_ParamName2 = param2_name_1,
        Outcome_1_Param2 = ifelse(length(est_1) > 1, round(est_1[2], 4), NA),
        Brunner_Munzel_P = round(pval, 5)
    )
    return(result)
}
variables <- c("lipids1", "lipids2", "lipids3", "lipids4", "lipids5", 
               "hormone1", "hormone2", "hormone3", "hormone4")
all_results <- do.call(rbind, lapply(variables, analyze_variable, data = data_for_analysis))
View(all_results)  # 或 print(all_results)
library(gtsummary)
tbl_summary(all_results, by="Variable")
library(DataExplorer)
create_report(all_results)  # Generates HTML report with graphs and statistics
create_report(
    data = all_results,
    output_file = "all_results.html",  
    output_dir = getwd(),                
    report_title = "all_results"          
)

# save data_for_analysis in CSV
write.csv(all_results, "all_results.csv", row.names = FALSE)  
