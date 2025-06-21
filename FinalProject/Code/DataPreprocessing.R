data <- read.csv("framingham.csv")

str(data)
summary(data)
sum(is.na(data))   
colSums(is.na(data))

# install.packages("skimr")
library(skimr)
skim(data)




# 安装可视化缺失值的包（如未安装）
install.packages("visdat")
install.packages("naniar")


na_stats <- colMeans(is.na(data)) * 100 # % missing data
na_stats
na_stats_filtered <- na_stats[na_stats <= 35] #  missing data <=35 %


library(visdat)
library(naniar)

# 可视化缺失模式
vis_miss(data)               # 带颜色展示缺失位置
gg_miss_var(data)            # 每个变量缺失的 bar 图

# ------------ 异常值 -----------------
# 加载库
library(dplyr)
library(ggplot2)
library(tidyr)
library(dbscan)

# ------------ Step 1：指定用于异常检测的变量 ------------
lof_variables <- c("BMI", "diaBP", "glucose", "heartRate")

# 选取数据并添加 row_id 用于追踪
lof_input <- data %>%
    select(all_of(lof_variables)) %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(cols = -row_id, names_to = "variable", values_to = "value")

heavy_smokers <- data %>%
    filter(cigsPerDay > 50)
heavy_glucose <- data %>%
    filter(glucose > 300)

# ------------ Step 2：绘制箱型图 ------------
ggplot(lof_input, aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    labs(title = "Outlier Detection (Boxplots)",
         x = "Variables",
         y = "Value") +
    theme_minimal()

# ------------ Step 3：使用 LOF 识别异常值 ------------
# 重新组织为原始格式（每行是样本）
lof_matrix <- data %>%
    select(all_of(lof_variables)) %>%
    na.omit()

# 保存行号（方便回溯）
row_ids <- which(complete.cases(data[, lof_variables]))

# 计算 LOF 分数
lof_scores <- lof(lof_matrix, k = 20)

# 构建带行号的结果框
lof_df <- data.frame(
    row_id = row_ids,
    lof_score = lof_scores
)

# 绘制直方图
ggplot(lof_df, aes(x = lof_score)) +
    geom_histogram(binwidth = 0.05, fill = "#FF7F00", color = "black", alpha = 0.7) +
    labs(title = "Histogram of LOF Scores", x = "LOF Score", y = "Frequency") +
    theme_minimal()

# ------------ Step 4：筛选异常值（top 5% 或 lof_score > 1.5） ------------
threshold <- quantile(lof_df$lof_score, 0.95)
lof_df <- lof_df %>%
    mutate(is_outlier = lof_score > threshold)

lof_suspects <- lof_df %>%
    filter(lof_score > 1.5) %>%
    arrange(desc(lof_score))

top_lof_errors <- lof_suspects %>%
    slice(1:50)

print(top_lof_errors)

# ------------ Step 5：z-score 验证异常值列 ------------
# 用于计算 z-score 的参考数据
z_reference <- data %>%
    select(all_of(lof_variables))

# 定义 z-score 检查函数
zscore_outlier_check <- function(row, ref_df, threshold = 2) {
    sapply(names(row), function(col) {
        mu <- mean(ref_df[[col]], na.rm = TRUE)
        sigma <- sd(ref_df[[col]], na.rm = TRUE)
        abs((row[[col]] - mu) / sigma) > threshold
    })
}

# 计算 top 异常行在各变量上的 z-score 异常
error_matrix <- t(sapply(top_lof_errors$row_id, function(i) {
    zscore_outlier_check(data[i, lof_variables], z_reference)
}))

# 构建结果矩阵
error_matrix_df <- as.data.frame(error_matrix)
rownames(error_matrix_df) <- paste0("Row_", top_lof_errors$row_id)

# 查看结果
View(error_matrix_df)


# ------------ 插值 ------------
# install.packages("mice")
library(mice)
library(dplyr)

# 对字符型变量转换为因子，避免报错
data_clean <- data %>%
    select(where(~!all(is.na(.)))) %>%
    mutate(across(where(is.character), as.factor))

# 执行 Little's MCAR 检验
mcar_result <- mcar_test(data_clean)
# 输出结果并解释
print(mcar_result)

# 自动解释函数
interpret_mcar <- function(mcar_result) {
    p <- mcar_result$p.value
    if (p > 0.05) {
        message("value > 0.05 → Data is likely MCAR. Safe to delete or impute.")
    } else {
        message("value <= 0.05 → Data is NOT MCAR. Assume MAR or MNAR.")
    }
}
interpret_mcar(mcar_result)

# imputed_data_rf <- mice(data, m=5, method='rf', print=FALSE)
imputed_data_rf <- mice(data[, !names(data) %in% "New"], method="rf")  
imputed_data_rf_final <- complete(imputed_data_rf)  # generate full data

library(ggplot2)
# Density plots 
ggplot(data, aes(x=glucose, fill="Original")) +
    geom_density(alpha=0.5) +
    geom_density(data=imputed_data_rf_final, aes(x=glucose, fill="Imputed"), alpha=0.5) +
    labs(title="Density Plot of glucose: Original vs. Imputed")

imputed_data_pmm <- mice(data[, !names(data) %in% "New"], method="pmm")  
imputed_data_pmm_final <- complete(imputed_data_pmm)  # generate full data
# Density plots 
ggplot(data, aes(x=glucose, fill="Original")) +
    geom_density(alpha=0.5) +
    geom_density(data=imputed_data_pmm_final, aes(x=glucose, fill="Imputed"), alpha=0.5) +
    labs(title="Density Plot of hormone10_generated: Original vs. Imputed")

# Camparsion
library(ggplot2)
library(dplyr)
library(visdat)

# 自动提取需要插值的变量名（有缺失的数值变量）
vars_to_plot <- data %>%
    select(where(is.numeric)) %>%
    select(where(~ any(is.na(.)))) %>%
    names()

# 定义统一的绘图函数
plot_imputation_comparison <- function(var_name) {
    # 提取三组数据
    original <- data %>%
        select(all_of(var_name)) %>%
        mutate(source = "Original")
    
    rf_imputed <- imputed_data_rf_final %>%
        select(all_of(var_name)) %>%
        mutate(source = "Random Forest")
    
    pmm_imputed <- imputed_data_pmm_final %>%
        select(all_of(var_name)) %>%
        mutate(source = "PMM")
    
    # 合并数据
    combined <- bind_rows(original, rf_imputed, pmm_imputed)
    
    if (var_name == "hormone10_generated") {
        x_lim <- c(0, 5)
    } else {
        # 自动获取范围
        x_min <- min(combined[[var_name]], na.rm = TRUE)
        x_max <- max(combined[[var_name]], na.rm = TRUE)
        
        # 加入buffer，防止贴边
        x_range <- x_max - x_min
        buffer <- x_range * 0.05
        x_lim <- c(x_min - buffer, x_max - buffer)
    }
    
    # 绘制密度图
    ggplot(combined, aes_string(x = var_name, fill = "source", color = "source")) +
        geom_density(alpha = 0.4, size = 1, na.rm = TRUE) +
        labs(title = paste("Density Comparison of:", var_name),
             x = var_name,
             y = "Density") +
        scale_x_continuous(limits = x_lim) +
        scale_fill_manual(values = c(
            "Original" = "black",
            "PMM" = "#E69F00",
            "Random Forest" = "#56B4E9"
        )) +
        scale_color_manual(values = c(
            "Original" = "black",
            "PMM" = "#E69F00",
            "Random Forest" = "#56B4E9"
        )) +
        theme_minimal() +
        theme(legend.position = "top")
}

for (v in vars_to_plot) {
    print(plot_imputation_comparison(v))     # 显示图形
    # 可选：保存图像（取消注释保存）
    ggsave(paste0("imputation_density_", v, ".png"), 
           plot = plot_imputation_comparison(v),
           width = 6, height = 4)
}
colSums(is.na(imputed_data_pmm_final))


smoking_logic_error <- data %>%
    filter(currentSmoker == 0 & cigsPerDay > 0)
print(smoking_logic_error)

write.csv(imputed_data_pmm_final, file = "imputed_data_pmm_final.csv", row.names = FALSE)