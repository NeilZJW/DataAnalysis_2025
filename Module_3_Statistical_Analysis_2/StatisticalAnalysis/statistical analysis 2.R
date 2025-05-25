#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset--------------------------
data <- read.csv("data_for_analysis_imputed.csv")
summary(data)
# testing for normality of distribution
shapiro.test(data$lipids1)
shapiro.test(data$lipids2)

hist(data$lipids1)  
qqnorm(data$lipids1)

# Spearman's correlation test

spearman_result<-cor.test(data$lipids1, data$lipids2, method="spearman")

print(spearman_result)

# data.frame for result
results <- data.frame(
  variable = character(),
  spearman_corr = numeric(),
  s_p_value = numeric(),
  stringsAsFactors = FALSE
)

# variables for analysis
target_vars <- c("lipids2", "lipids3", "lipids4")

library(wPerm)

# main 
for (var in target_vars) {
  # spearman
  perm_spearman <- perm.relation(
    x = data$lipids1, 
    y = data[[var]],
    method = "spearman",
    R = 10000
  )
  
  # add result
  results <- rbind(results, data.frame(
    variable = var,
    spearman_corr = perm_spearman$Observed,
    s_p_value =  perm_spearman$p.value
    ))
}


# output result
print(results)

#------visualization of significant results of correlation analysis---------

data<-data[order(data$lipids1),]

plot(data$lipids1, data$lipids2)

lines(data$lipids1, data$lipids2, col = "blue")

abline(lm(data$lipids1 ~ data$lipids2), col="red")




#_____________regression analysis________________ 

df=data
df<-df[order(df$lipids1),]


#linear regression

model_linear <- lm(lipids1 ~ lipids2, data=df)
summary(model_linear)


#second degree polynomal

model_2 <- lm(lipids1 ~ poly(lipids2, 2), data=df)
summary(model_2)

#third degree polynomal

model_3 <- lm(lipids1 ~ poly(lipids2, 3), data=df)

summary(model_3)
#exponential dependence

model_exp <- lm(log(lipids1) ~ lipids2, data=df)
summary(model_exp)
# log dependence

model_log <- lm(exp(lipids1) ~ lipids2, data=df)
summary(model_log)
#comparison of models
#table of result

rezult<-data.frame(model=c("model_linear", "model_2", "model_3", "model_exp", "model_log"), BIC_value=c(BIC(model_linear), BIC(model_2), BIC(model_3), BIC(model_exp), BIC(model_log)))

rezult<-rezult[order(rezult$BIC_value),]

rezult


# __________building graphs______________
#         linear regression graphs

plot(df$lipids1, df$lipids2)


# HomeWork

data <- read.csv("data_for_analysis_imputed.csv")
numeric_vars <- sapply(data, is.numeric)
numeric_data <- data[, numeric_vars]

target <- "lipids5"

other_vars <- setdiff(names(numeric_data), target)


# 创建结果表格
shapiro_results <- data.frame(
    variable = character(),
    p_value = numeric(),
    is_normal = logical(),
    stringsAsFactors = FALSE
)

# 对每个数值型变量进行 Shapiro-Wilk 检验
for (var in names(numeric_data)) {
    values <- numeric_data[[var]]
    
    # 去除缺失值
    values <- values[!is.na(values)]
    
    # Shapiro-Wilk 检验（样本数必须 < 5000）
    if (length(values) >= 3 && length(values) < 5000) {
        test <- shapiro.test(values)
        shapiro_results <- rbind(shapiro_results, data.frame(
            variable = var,
            p_value = test$p.value,
            is_normal = test$p.value >= 0.05
        ))
    } else {
        shapiro_results <- rbind(shapiro_results, data.frame(
            variable = var,
            p_value = NA,
            is_normal = NA
        ))
    }
}

# 打印结果
print(shapiro_results)
write.csv(shapiro_results, "shapiro_wiki_results.csv", row.names = FALSE)



na_columns <- names(data)[colSums(is.na(data)) > 0]
data_clean <- data[, !(names(data) %in% na_columns)]
numeric_data <- data_clean[sapply(data_clean, is.numeric)]
print(names(numeric_data))
write.csv(numeric_data, "data_clean.csv", row.names = FALSE)


results <- data.frame(
    variable = character(),
    spearman_corr = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
)
for (var in other_vars) {
    tryCatch({
        # spearman
        perm_spearman <- perm.relation(
            x = numeric_data[[target]], 
            y = numeric_data[[var]],
            method = "spearman",
            R = 10000
        )
        
        # add result
        results <- rbind(results, data.frame(
            variable = var,
            spearman_corr = perm_spearman$Observed,
            s_p_value =  perm_spearman$p.value
        ))
    }, error = function(e) {
        message(paste("Error：", var))
    })
}
results <- results[order(results$s_p_value), ]
print(results)
write.csv(results, "spearman.csv", row.names = FALSE)

names(df)
# Attention: "сarb_metabolism" c is not English letter, that is Russian с
selected_vars <- c("lipids1", "lipids2", "lipids3", "сarb_metabolism", "hormone2", "lipid_pero4")
df <- read.csv("data_for_analysis_imputed.csv")

library(dplyr)
final_result <- data.frame()

for (var in selected_vars) {
    df_sorted <- df[order(df[[var]]), ]
    
    # 线性
    m1 <- lm(lipids5 ~ df_sorted[[var]], data=df_sorted)
    bic1 <- BIC(m1)
    
    # 二次多项式
    m2 <- lm(lipids5 ~ poly(df_sorted[[var]], 2), data=df_sorted)
    bic2 <- BIC(m2)
    
    # 三次多项式
    m3 <- lm(lipids5 ~ poly(df_sorted[[var]], 3), data=df_sorted)
    bic3 <- BIC(m3)
    
    # 指数回归：log(Y) ~ X
    if (all(df_sorted$lipids5 > 0)) {
        m4 <- lm(log(lipids5) ~ df_sorted[[var]], data=df_sorted)
        bic4 <- BIC(m4)
    } else {
        bic4 <- NA
    }
    
    # 对数回归：Y ~ log(X)
    if (all(df_sorted[[var]] > 0)) {
        m5 <- lm(lipids5 ~ log(df_sorted[[var]]), data=df_sorted)
        bic5 <- BIC(m5)
    } else {
        bic5 <- NA
    }
    
    temp <- data.frame(
        variable = var,
        model = c("linear", "poly2", "poly3", "exp", "log"),
        BIC = c(bic1, bic2, bic3, bic4, bic5)
    )
    
    final_result <- bind_rows(final_result, temp)
}

best_models <- final_result %>%
    group_by(variable) %>%
    slice_min(BIC, n = 1, with_ties = FALSE)

print(final_result)
print(best_models)
write.csv(final_result, "regression_BIC_all_models.csv", row.names = FALSE)
write.csv(best_models, "regression_BIC_best_models.csv", row.names = FALSE)

for (var in selected_vars) {
    data_sorted <- df[order(df[[var]]), ]
    
    # 设置保存文件路径
    file_name <- paste0("plots/lipids5_vs_", var, ".png")
    
    # 打开图像保存设备
    png(filename = file_name, width = 800, height = 600)
    
    # 绘图 + 拟合线（示例风格）
    plot(data_sorted[[var]], data_sorted$lipids5,
         main = paste("lipids5 vs", var),
         xlab = var, ylab = "lipids5",
         pch = 16, col = "black")
    
    lines(data_sorted[[var]], data_sorted$lipids5, col = "blue", lwd = 1)
    
    abline(lm(data_sorted$lipids5 ~ data_sorted[[var]]), col = "red", lwd = 2)
    
    # 关闭图像保存设备
    dev.off()
}

