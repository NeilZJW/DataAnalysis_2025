data <- read.csv("imputed_data_pmm_final.csv")

str(data)
summary(data)
sum(is.na(data))   
colSums(is.na(data))

# testing for normality of distribution
shapiro.test(data$TenYearCHD)
shapiro.test(data$glucose)

hist(data$glucose)  
qqnorm(data$glucose)


# Spearman's correlation test

spearman_result<-cor.test(data$TenYearCHD, data$glucose, method="spearman")
print(spearman_result)

# 加载包
library(wPerm)

# 提取数值型变量名（去除 TenYearCHD 自身）
numeric_vars <- names(data)[sapply(data, is.numeric)]
target_vars <- setdiff(numeric_vars, "TenYearCHD")

# 结果存储表
spearman_results <- data.frame(
    variable = character(),
    spearman_corr = numeric(),
    s_p_value = numeric(),
    stringsAsFactors = FALSE
)

# 遍历每个变量与 TenYearCHD 做 Spearman 检验
for (var in target_vars) {
    test_result <- perm.relation(
        x = data[[var]],
        y = data$TenYearCHD,
        method = "spearman",
        R = 10000
    )
    
    spearman_results <- rbind(spearman_results, data.frame(
        variable = var,
        spearman_corr = test_result$Observed,
        s_p_value = test_result$p.value
    ))
}

# 按照相关性大小排序
spearman_results <- spearman_results[order(abs(spearman_results$spearman_corr), decreasing = TRUE), ]

# 打印前几个结果
print(spearman_results)

# 可选：保存结果为 CSV
write.csv(spearman_results, "spearman_results.csv", row.names = FALSE)


#------visualization of significant results of correlation analysis---------
library(ggplot2)

ggplot(data, aes(x = factor(TenYearCHD), y = age, fill = factor(TenYearCHD))) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "darkred") +
    labs(title = "Age vs. CHD", x = "TenYearCHD", y = "Age") +
    theme_minimal()


# --------------- 逻辑回归 ----------------
# 加载必要的包
library(car)

# 构建逻辑回归模型（family = binomial 表示逻辑回归）
model_logistic <- glm(
    TenYearCHD ~ age + sysBP + prevalentHyp + diaBP + diabetes + BPMeds + male + totChol + BMI + education + prevalentStroke + glucose + cigsPerDay,
    data = data,
    family = binomial
)

# 查看模型摘要
summary(model_logistic)
model_logistic$converged
vif(model_logistic)

step_model <- step(model_logistic, direction = "both")
summary(step_model)



install.packages("pROC")

# 然后加载
library(pROC)


pred_old <- predict(model_logistic, type = "response")
pred_new <- predict(step_model, type = "response")

# 绘图对比
roc_old <- roc(data$TenYearCHD, pred_old)
roc_new <- roc(data$TenYearCHD, pred_new)

plot(roc_old, col = "red")
lines(roc_new, col = "blue")
legend("bottomright", legend = c("Full Model", "Stepwise Model"), col = c("red", "blue"), lwd = 2)

auc(roc_old)
auc(roc_new)
table(data$TenYearCHD, pred_old > 0.5)
table(data$TenYearCHD, pred_new > 0.5)

install.packages("sjPlot")
library(sjPlot)
plot_model(step_model, type = "est", show.values = TRUE)

table(data$TenYearCHD, pred_new > 0.4)

