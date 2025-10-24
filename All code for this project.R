library(tidyverse)#数据处理包
library(haven)#加载dta格式数据
mimic_data <- read_dta("MIMIC-IV.dta")
eicu_data <- read_dta("eICU-CRD.dta")

#观察缺失值
library(naniar)
missing_summary_mimic <- miss_var_summary(mimic_data)#显示各列缺失值占比
print(missing_summary_mimic)
gg_miss_var(mimic_data) + theme_bw()#缺失值可视化
missing_summary_eicu <- miss_var_summary(eicu_data)#显示各列缺失值占比
print(missing_summary_eicu)
gg_miss_var(eicu_data) + theme_bw()#缺失值可视化

#生成生存分析及住院日箱线图数据未发现缺失值

#生成患者住院 30天院内死亡率生存时间数据（KM图）
library(survival)
library(survminer)
mimic_km <- mimic_data %>% select(sad,icustay,icu28dmort,hosp_mort,hospstay)
mimic_km <- mimic_km %>% mutate(icustay_d = round(icustay/60))
mimic_km <- mimic_km %>% mutate(hospstay_d = round(hospstay/60))
mimic_km <- mimic_km %>% mutate(hospstay_30d = ifelse(hospstay_d >= 30,30,hospstay_d))
mimic_km <- mimic_km %>% mutate(hosp_30dmort = case_when(
  hospstay_d > 30 & hosp_mort == 1 ~ 0,
  hospstay_d > 30 & hosp_mort == 0 ~ 0,
  hospstay_d <= 30 & hosp_mort == 1 ~ 1,
  hospstay_d <= 30 & hosp_mort == 0 ~ 0
))#建立新的30天院内死亡率数据

diff <- survdiff(Surv(hospstay_30d,hosp_30dmort)~sad,data = mimic_km)
pvalue <- 1-pchisq(diff$chisq,df=(2-1))#df自由度
if(pvalue<0.001){
  pvalue="Log-rank test \n p<0.001"
}else{
  pvalue=paste0("Log-rank test","\n","p=",sprintf("%.03f",pvalue))
}
fit <- survfit(Surv(hospstay_30d,hosp_30dmort)~sad,data = mimic_km)
#绘制
mimic_hosp_30d <- ggsurvplot(fit,
                             data = mimic_km,
                             conf.int = F,
                             #pval = T,
                             pval = pvalue,
                             pval.size=5,
                             pval.coord = c(0, 0.5),  # 调整P值位置，x=0, y=0.5
                             legend.labs=c("Non-SAD","SAD"),
                             #legend.labs=levels(factor(mimic_km[,"sad"])),
                             legend.title="Group",
                             xlab="Time(day)",
                             break.time.by=5,
                             risk.table.title="Number at risk",
                             risk.table = T,
                             risk.table.height=.25,
                             censor = T,
                             censor.shape=124,
                             censor.size=2,
                             ggtheme = theme_bw(),
                             ylim = c(0.25, 1),
                             type = "spline",  # 使生存曲线更圆滑
                             palette = c("blue", "red")  # 设置曲线颜色
                             )  
print(mimic_hosp_30d)
jpeg("Figure 2A survival_plot.jpeg", width = 8, height = 6, units = "in", res = 600)
print(mimic_hosp_30d)
dev.off()

#生存住院日箱线图
library(tidyverse)
library(ggpubr)
mimic_boxplot <- mimic_data %>% select(sad,hosp_mort,hospstay)
mimic_boxplot <- mimic_boxplot %>% mutate(hospstay_d = round(hospstay/60))
mimic_boxplot <- mimic_boxplot %>% mutate(hospstay_90d = ifelse(hospstay_d >=50,50,hospstay_d))

mimic_boxplot$sad <- factor(mimic_boxplot$sad,
                                    levels = c("0","1"),
                                    labels = c("Non-SAD","SAD"))#将分组设置标签
mimic_hosp_d <- ggboxplot(mimic_boxplot,x="sad",y="hospstay_90d",color = "sad",
          palette = c("blue","red"),
          ggtheme = theme_bw(),
          legend.title="Group",
          legend.labs=c("Non-SAD","SAD"),
          xlab = "",
          ylab = "Time(day)",
          ylim = c(0,50)
          #main= "MIMIC-IV 数据库中谵妄组和非谵妄组的 ICU 住院时间箱线图"
)+
  #stat_compare_means()#添加Wilcoxon秩和检验计算的P值+
  #stat_compare_means(label = "p.signif",
  #symnum.args = list(cutpoints=c(0,0.0001,0.001,0.01,0.05,Inf),
  #symbols=c("<0.0001","<0.001","<0.01","<0.05","ns")))
  stat_compare_means(label = "p.signif",
                     symnum.args = list(cutpoints=c(0,0.01,0.05,Inf),
                                        symbols=c("P<0.01","P<0.05","ns")),#获取区间并设置区间的P值怎么显示
                     label.x.npc = "right",#确定P值的x轴位置
                     label.y.npc = "top")#确定P值的y轴位置
print(mimic_hosp_d)
jpeg("Figure 2B boxplot.jpeg", width = 8, height = 6, units = "in", res = 600)
print(mimic_hosp_d)
dev.off()



#患者30天院内死亡率生存组与死亡组基线对比
library(tidyverse)
#整理数据并对缺失值进行多重插补

#筛选SAD患者
mimic_sad <- mimic_data %>% filter(sad==1)
eicu_sad <- eicu_data %>% filter(sad==1)

#对缺失值进行数据插补
library(mice)#使用mice包进行数据多重插补
set.seed(123)
mimic_sad <- mimic_sad[,-1]#多重插补不能有患者ID以及编号之类的数据，需要剔除
mimic_sad <- mice(mimic_sad,m=5,maxit = 10,method = "pmm",seed = 123,print = TRUE)#设置种子，输出过程
mimic_sad <- mice::complete(mimic_sad)
eicu_sad <- eicu_sad[,-1]#多重插补不能有患者ID以及编号之类的数据，需要剔除
eicu_sad <- mice(eicu_sad,m=5,maxit = 10,method = "pmm",seed = 123,print = TRUE)#设置种子，输出过程
eicu_sad <- mice::complete(eicu_sad)

#整理表格内住院30天院内死亡数据
#建立MIMIC-IV新的30天院内死亡率数据
mimic_sad <- mimic_sad %>% mutate(hospstay_d = round(hospstay/60))
mimic_sad <- mimic_sad %>% mutate(hospstay_30d = ifelse(hospstay_d >= 30,30,hospstay_d))
mimic_sad <- mimic_sad %>% mutate(hosp_30dmort = case_when(
  hospstay_d > 30 & hosp_mort == 1 ~ 0,
  hospstay_d > 30 & hosp_mort == 0 ~ 0,
  hospstay_d <= 30 & hosp_mort == 1 ~ 1,
  hospstay_d <= 30 & hosp_mort == 0 ~ 0
))

#建立eICU-CRD新的30天院内死亡率数据
eicu_sad <- eicu_sad %>% mutate(hospstay_d = round(hospstay/60))
eicu_sad <- eicu_sad %>% mutate(hospstay_30d = ifelse(hospstay_d >= 30,30,hospstay_d))
eicu_sad <- eicu_sad %>% mutate(hosp_30dmort = case_when(
  hospstay_d > 30 & hosp_mort == 1 ~ 0,
  hospstay_d > 30 & hosp_mort == 0 ~ 0,
  hospstay_d <= 30 & hosp_mort == 1 ~ 1,
  hospstay_d <= 30 & hosp_mort == 0 ~ 0
))


#整理数据（以hosp_30dmort为应变量（放在第一列），减少不需要的临床特征）
mimic_sad <- mimic_sad %>% select(-c(race,first_careunit,deliriumtime,sepsistime,icu28dmort)) %>% select(hosp_30dmort,everything())
eicu_sad <- eicu_sad %>% select(-c(race,first_careunit,deliriumtime,sepsistime,icu28dmort)) %>% select(hosp_30dmort,everything())

#生成表1
library(tableone)#生成基线对比表格
#分类变量转换
dput(names(mimic_sad))#获取变量名
#纳入变量定义（排除ID及不需要的变量）
myvars <- c("hosp_30dmort", "age", "weight", "gender", "temperature", "heart_rate", 
            "resp_rate", "spo2", "sbp", "dbp", "mbp", "wbc", "hemoglobin", 
            "platelet", "bun", "cr", "glu", "Na", "Cl", "K", "Mg", "Ca", 
            "P", "inr", "pt", "ptt", "bicarbonate", "aniongap", "gcs", "vent", 
            "crrt", "vaso", "seda", "sofa_score", "ami", "ckd", "copd", "hyperte", 
            "dm", "aki", "stroke")
#分类变量定义
catvars <- c("hosp_30dmort", "gender","vent", 
             "crrt", "vaso", "seda","ami", "ckd", "copd", "hyperte", 
             "dm", "aki", "stroke")
#分组汇总+偏态数据使用中位数表示
#生成MIMIC-IV数据库表1
tab_mimic <- CreateTableOne(vars = myvars,strata = "hosp_30dmort",data = mimic_sad,factorVars = catvars)
tableone_mimic <- print(tab_mimic,nonnormal = TRUE,quote = FALSE,noSpaces = TRUE,printToggle = FALSE)
write.csv(tableone_mimic,file="tableone_mimic.csv")
#生成eICU-CRD数据框表1
tab_eicu <- CreateTableOne(vars = myvars,strata = "hosp_30dmort",data = eicu_sad,factorVars = catvars)
tableone_eicu <- print(tab_eicu,nonnormal = TRUE,quote = FALSE,noSpaces = TRUE,printToggle = FALSE)
write.csv(tableone_eicu,file="tableone_eicu.csv")

#对所有变量进行单变量逻辑回归
mimic_sad_LR <- mimic_sad %>% select(hosp_30dmort, age, weight, gender, temperature, heart_rate, 
                                     resp_rate, spo2, sbp, dbp, mbp, wbc, hemoglobin, 
                                     platelet, bun, cr, glu, Na, Cl, K, Mg, Ca, 
                                     P, inr, pt, ptt, bicarbonate, aniongap, gcs, vent, 
                                     crrt, vaso, seda, sofa_score, ami, ckd, copd, hyperte, 
                                     dm, aki, stroke)
# 初始化结果列表
results_mimic_LR_list <- list()

# 遍历每个自变量进行单因素逻辑回归
for (i in 2:ncol(mimic_sad_LR)) {
  # 构建模型
  model_LR <- glm(hosp_30dmort ~ mimic_sad_LR[, i], data = mimic_sad_LR, family = binomial)
  # 存储结果
  results_mimic_LR_list[[colnames(mimic_sad_LR)[i]]] <- summary(model_LR)
}

# 查看结果
for (var_name in names(results_mimic_LR_list)) {
  cat("变量:", var_name, "\n")
  print(results_mimic_LR_list[[var_name]]$coefficients)
  cat("\n")
}

#单变量逻辑回归OR值及置信区间提取
mimic_sad_LR <- mimic_sad %>% select(hosp_30dmort, age, weight, gender, temperature, heart_rate, 
                                     resp_rate, spo2, sbp, dbp, mbp, wbc, hemoglobin, 
                                     platelet, bun, cr, glu, Na, Cl, K, Mg, Ca, 
                                     P, inr, pt, ptt, bicarbonate, aniongap, gcs, vent, 
                                     crrt, vaso, seda, sofa_score, ami, ckd, copd, hyperte, 
                                     dm, aki, stroke)
# 初始化结果列表
results <- list()

# 遍历每个自变量进行单因素逻辑回归
for (i in 2:ncol(mimic_sad_LR)) {
  # 构建模型
  model <- glm(hosp_30dmort ~ mimic_sad_LR[, i], data = mimic_sad_LR, family = binomial)
  
  # 提取回归系数和置信区间
  coefficients <- summary(model)$coefficients
  ci <- exp(confint(model))  # 置信区间指数运算
  
  # 存储结果
  results[[colnames(mimic_sad_LR)[i]]] <- data.frame(
    OR = exp(coefficients[, 1]),  # OR 值
    Lower_CI = ci[, 1],           # 置信区间下限
    Upper_CI = ci[, 2],           # 置信区间上限
    P_value = coefficients[, 4]   # p 值
  )
}

# 查看结果
for (var_name in names(results)) {
  cat("变量:", var_name, "\n")
  print(results[[var_name]])
  cat("\n")
}




#进行相关性分析以减少多重共线性
#生成相关性分析数据框
library(corrplot)
mimic_sad_spearman <- mimic_sad_LR %>% select(age, weight, gender, temperature, heart_rate, 
                                              resp_rate, spo2, sbp, dbp, mbp, wbc, hemoglobin, 
                                              bun, cr, Cl, K, Mg, 
                                              P, inr, pt, ptt, bicarbonate, aniongap,
                                              crrt, vaso, seda, sofa_score, ami, ckd, hyperte, 
                                              aki, stroke)
# 计算 Spearman 相关系数矩阵
spearman_cor_matrix <- cor(mimic_sad_spearman, method = "spearman", use = "complete.obs")
# 查看相关系数矩阵
print(spearman_cor_matrix)

# 可视化相关系数矩阵（可选）
jpeg("Figure 6 Spearman heatmap plot.jpeg", width = 8, height = 6, units = "in", res = 600)
corrplot(spearman_cor_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

#筛选相关系数绝对值大于0.75的变量对
significant_cor <- which(abs(spearman_cor_matrix) > 0.75, arr.ind = TRUE)
# 查看结果
print(significant_cor)
#导出相关系数
write.csv(as.data.frame(spearman_cor_matrix),file="mimic_sad_spearman.csv")



#根据单因素回归结果进行Boruta算法筛选特征值
library(Boruta)
library(tidyverse)
mimic_sad_boruta <- mimic_sad_LR %>% select(hosp_30dmort, age, weight, gender, temperature, heart_rate, 
                                        resp_rate, spo2, mbp, wbc, hemoglobin, 
                                        bun, Cl, K, Mg, 
                                        P, inr, ptt, bicarbonate, aniongap,
                                        crrt, vaso, seda, sofa_score, ami, ckd, hyperte, 
                                        aki, stroke)
str(mimic_sad_boruta)
mimic_sad_boruta$hosp_30dmort <- as.factor(mimic_sad_boruta$hosp_30dmort)

set.seed(123)  # 设置随机种子以确保结果可重复
boruta_output <- Boruta(hosp_30dmort ~ ., data = mimic_sad_boruta, doTrace = 2, maxRuns = 300)
print(boruta_output)
plot(boruta_output, las = 2, cex.axis = 0.8)#绘制结果
jpeg("Figure 3 Boruta.jpeg", width = 8, height = 6, units = "in", res = 600)
plot(boruta_output, las = 2, cex.axis = 0.6)#绘制结果
dev.off()
plotImpHistory(boruta_output)#查看特征重要性历史
getSelectedAttributes(boruta_output)#获取确认的特征名称
getConfirmedFormula(boruta_output)#获取确认的特征的公式
selected_features <- getSelectedAttributes(boruta_output)

#根据Boruta算法结果筛选特征变量
mimic_ml <- mimic_sad_boruta %>% dplyr::select(hosp_30dmort,all_of(selected_features))


missing_summary_mimic_ml <- miss_var_summary(mimic_ml)
print(missing_summary_mimic_ml)#观察是否存在缺失值

library(caret)
mimic_sad_ml <- mimic_ml
set.seed(123)
index <- createDataPartition(mimic_sad_ml$hosp_30dmort,p=0.7,list = F)
train <- mimic_sad_ml[index,]
test <- mimic_sad_ml[-index,]

#将训练集中的相应列转换为因子类型
train$hosp_30dmort <- as.factor(as.character(train$hosp_30dmort))
train$vaso <- as.factor(train$vaso)
train$seda <- as.factor(train$seda)
train$crrt <- as.factor(train$crrt)
train$ckd <- as.factor(train$ckd)
train$aki <- as.factor(train$aki)
train$stroke <- as.factor(train$stroke)
str(train)

#将内部验证集中的相应列转换为因子类型
test$hosp_30dmort <- as.factor(as.character(test$hosp_30dmort))
test$vaso <- as.factor(test$vaso)
test$seda <- as.factor(test$seda)
test$crrt <- as.factor(test$crrt)
test$ckd <- as.factor(test$ckd)
test$aki <- as.factor(test$aki)
test$stroke <- as.factor(test$stroke)
str(test)

#LR model
lm_model <- glm(hosp_30dmort~.,data = train,family = binomial(link = "logit"))
summary(lm_model)
lm_pred <- predict(lm_model,test,type="response")
threshold <- 0.5
predictions_binary <- ifelse(lm_pred>threshold,1,0)
confusionMatrix(as.factor(predictions_binary),as.factor(test$hosp_30dmort))
LR_pred <- predictions_binary
summary(lm_model)#获取模型的详细信息
or_values <- exp(coef(lm_model))# 计算 OR 值
ci <- confint(lm_model)# 计算系数的置信区间
or_ci <- exp(ci)# 计算 OR 值的置信区间
# 组合结果并输出
or_results <- data.frame(
  OR = or_values,
  lower_CI = or_ci[, 1],
  upper_CI = or_ci[, 2]
)
print(or_results)
#lm_f1 <- confusionMatrix(as.factor(predictions_binary),as.factor(test$hosp_30dmort))$byClass[2]

#SVM model
library(e1071)
svm_model <- svm(hosp_30dmort ~ ., data = train, probability = TRUE)
svm_pred <- predict(svm_model, test, probability = TRUE)
a <- data.frame(svm_pred)
ab <- as.factor(as.character(a$svm_pred))
svm_hosp_30dmort <- as.factor(test$hosp_30dmort)
svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]
confusionMatrix(data = ab, reference = svm_hosp_30dmort)
train$hosp_30dmort <- as.numeric(as.character(train$hosp_30dmort)) 
test$hosp_30dmort <- as.numeric(as.character(test$hosp_30dmort)) 
SVM_pred <- as.numeric(a$svm_pred)
svm_f1 <- confusionMatrix(data = ab, reference = svm_hosp_30dmort)$byClass[2]

#XGBoost model
train$hosp_30dmort <- as.numeric(as.character(train$hosp_30dmort))#确保应变量为数值类型
test$hosp_30dmort <- as.numeric(as.character(test$hosp_30dmort))#确保应变量为数值类型
library(xgboost)
train_matrix <- xgb.DMatrix(data.matrix(train[,-1]),label = train$hosp_30dmort)
test_matrix <- xgb.DMatrix(data.matrix(test[,-1]),label = test$hosp_30dmort)
params <- list(objective = "binary:logistic",eval_metric="logloss",max_depth=3,eta = 0.1,
               gamma=0.5,colsample_bytree = 1,min_child_weight=1,subsample=0.5)
watchlist <- list(train = train_matrix,val=test_matrix)
xgb_model <- xgb.train(params = params,data = train_matrix,nrounds = 125,watchlist = watchlist,
                       early_stopping_rounds = 10,print_every_n = 10,maximize = FALSE)
xgb_pred_prob <- predict(xgb_model,test_matrix)
xgb_pred <- ifelse(xgb_pred_prob > 0.5,1,0)
xgb_pred_factor <- factor(xgb_pred,levels = c(0,1))
test_hosp_30dmort_factor <- factor(test$hosp_30dmort,levels = c(0,1))
confusionMatrix(data = xgb_pred_factor,reference=test_hosp_30dmort_factor)
#xgb_f1 <- confusionMatrix(data = xgb_pred_factor,reference=test_hosp_30dmort_factor)$byClass[2]


#RF model（随机森林）
library(randomForest)
set.seed(123)
train$hosp_30dmort <- as.factor(train$hosp_30dmort)
test$hosp_30dmort <- as.factor(test$hosp_30dmort)
rf_model <- randomForest(hosp_30dmort~.,data = train,ntree = 500,mtry = 6)
rf_pred <- predict(rf_model,newdata = test)
confusionMatrix(data = rf_pred,reference = test$hosp_30dmort)
#rf_f1 <- confusionMatrix(data = rf_pred,reference = test$hosp_30dmort)$byClass[2]

#DT model(决策树)
library(rpart)
dt_model <- rpart(hosp_30dmort~.,data = train,method = "class")
dt_pred_prob <- predict(dt_model,newdata = test,type="prob")[,2]
dt_pred <- ifelse(dt_pred_prob > 0.5,1,0)
confusionMatrix(factor(dt_pred,levels = c("0","1")),test$hosp_30dmort)
#dt_f1 <- confusionMatrix(factor(dt_pred,levels = c("0","1")),test$hosp_30dmort)$byClass[2]

#NB model(朴素贝叶斯)
library(e1071)
nb_model <- naiveBayes(hosp_30dmort~.,data = train)
nb_pred_prob <- predict(nb_model,newdata = test,type="raw")[,2]
nb_pred <- ifelse(nb_pred_prob > 0.5,1,0)
confusionMatrix(factor(nb_pred,levels = c("0","1")),test$hosp_30dmort)
#nb_f1 <- confusionMatrix(factor(nb_pred,levels = c("0","1")),test$hosp_30dmort)$byClass[2]

#KNN model(K最近邻回归)
library(kknn)
knn_model <- kknn(hosp_30dmort~.,train,test,k=10,distance = 2,kernel = "rectangular")
knn_pred_prob <- predict(knn_model,newdata=test,type="prob")
knn_pred_prob <- knn_pred_prob[,"1"]
knn_pred_prob <- as.numeric(knn_pred_prob)
threshold <- 0.5
knn_pred <- ifelse(knn_pred_prob > threshold,1,0)
confusionMatrix(factor(knn_pred,levels = c("0","1")),test$hosp_30dmort)
#knn_f1 <- confusionMatrix(factor(knn_pred,levels = c("0","1")),test$hosp_30dmort)$byClass[2]

#内部验证ROC图绘制(pROC)
library(pROC)
library(ggplot2)
test_hosp_30dmort <- test$hosp_30dmort
ML_ROC <- data.frame(test_hosp_30dmort,
                     lm_pred,LR_pred,
                     svm_pred,svm_pred_prob,
                     xgb_pred_prob,xgb_pred,
                     rf_pred,
                     dt_pred_prob,dt_pred,
                     nb_pred_prob,nb_pred,
                     knn_pred,knn_pred_prob
)
ML_ROC$test_hosp_30dmort <- as.numeric(as.character(ML_ROC$test_hosp_30dmort))
ML_ROC$svm_pred <- as.numeric(as.character(ML_ROC$svm_pred))
ML_ROC$rf_pred <- as.numeric(as.character(ML_ROC$rf_pred))
str(ML_ROC)
#指定结局变量并命名
LR_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$lm_pred);LR_ROC
SVM_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$svm_pred_prob);SVM_ROC
XGB_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$xgb_pred_prob);XGB_ROC
RF_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$rf_pred);RF_ROC
DT_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$dt_pred_prob);DT_ROC
NB_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$nb_pred_prob);NB_ROC
KNN_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$knn_pred_prob);KNN_ROC

#输出AUC的95%置信区间
auc(LR_ROC)
auc(SVM_ROC)
auc(XGB_ROC)
auc(RF_ROC)
auc(DT_ROC)
auc(NB_ROC)
auc(KNN_ROC)
ci.auc(LR_ROC)
ci.auc(SVM_ROC)
ci.auc(XGB_ROC)
ci.auc(RF_ROC)
ci.auc(DT_ROC)
ci.auc(NB_ROC)
ci.auc(KNN_ROC)

#输出每个模型的F1评分
#cat("LR 模型 F1 评分：", lm_f1, "\n")
#cat("SVM 模型 F1 评分：", svm_f1, "\n")
#cat("XGBoost 模型 F1 评分：", xgb_f1, "\n")
#cat("RF 模型 F1 评分：", rf_f1, "\n")
#cat("DT 模型 F1 评分：", dt_f1, "\n")
#cat("NB 模型 F1 评分：", nb_f1, "\n")
#cat("KNN 模型 F1 评分：", knn_f1, "\n")




#指定图像保存路径和文件名
#png("ML ROC curves.png",width = 800,height = 800)
jpeg(filename = "Figure 4A 内部验证ROC.jpeg", width = 6, height = 6, units = "in", res = 600)
#画出第一条曲线
plot.roc(SVM_ROC,
         max.auc.polygon=F,#填充整个图像
         smooth=F,#绘制不平滑曲线“
         #main="Comparison of different machine learning models of ROC curves",#添加标题
         col = "red",#曲线颜色为红色
         legacy.axes=T,#使横轴从0到1，表示为1-特异度
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.8
)
#逐步添加其他曲线
#plot.roc(SVM_ROC,
#         add = T,
#         col = "orange",
#         smooth = F,
#         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.75
#)
plot.roc(XGB_ROC,
         add = T,
         col = "black",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.70
)
plot.roc(RF_ROC,
         add = T,
         col = "pink",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.65
)
plot.roc(DT_ROC,
         add = T,
         col = "green",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.60
)
plot.roc(NB_ROC,
         add = T,
         col = "blue",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.55
)
plot.roc(KNN_ROC,
         add = T,
         col = "orange",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.50
)
#增加图例
legend("bottomright",
       legend = c("SVM（AUC=0.6996）","XGBoost（AUC=0.7308）","RF（AUC=0.5295）","DT（AUC=0.6388）","NB（AUC=0.6873）","KNN（AUC=0.6645）"),
       col = c("red","black","pink","green","blue","orange"),
       lwd = 1,
       cex=0.7
)
dev.off()

#外部验证ROC图绘制(pROC)
#根据Boruta算法结果筛选特征变量
eicu_ml <- eicu_sad %>% dplyr::select(hosp_30dmort,all_of(selected_features))


missing_summary_eicu_ml <- miss_var_summary(eicu_ml)
print(missing_summary_eicu_ml)#观察是否存在缺失值

library(caret)
set.seed(123)
test <- eicu_ml
#将训练集中的相应列转换为因子类型
train$hosp_30dmort <- as.factor(as.character(train$hosp_30dmort))
train$vaso <- as.factor(train$vaso)
train$seda <- as.factor(train$seda)
train$crrt <- as.factor(train$crrt)
train$ckd <- as.factor(train$ckd)
train$aki <- as.factor(train$aki)
train$stroke <- as.factor(train$stroke)
str(train)

#将外部验证集中的相应列转换为因子类型
test$hosp_30dmort <- as.factor(as.character(test$hosp_30dmort))
test$vaso <- as.factor(test$vaso)
test$seda <- as.factor(test$seda)
test$crrt <- as.factor(test$crrt)
test$ckd <- as.factor(test$ckd)
test$aki <- as.factor(test$aki)
test$stroke <- as.factor(test$stroke)
str(test)

#LR model
lm_model <- glm(hosp_30dmort~.,data = train,family = binomial(link = "logit"))
summary(lm_model)
lm_pred <- predict(lm_model,test,type="response")
threshold <- 0.5
predictions_binary <- ifelse(lm_pred>threshold,1,0)
confusionMatrix(as.factor(predictions_binary),as.factor(test$hosp_30dmort))
LR_pred <- predictions_binary


#SVM model
library(e1071)
svm_model <- svm(hosp_30dmort ~ ., data = train, probability = TRUE)
svm_pred <- predict(svm_model, test, probability = TRUE)
a <- data.frame(svm_pred)
ab <- as.factor(as.character(a$svm_pred))
svm_hosp_30dmort <- as.factor(test$hosp_30dmort)
svm_pred_prob <- attr(svm_pred, "probabilities")[, 2]
confusionMatrix(data = ab, reference = svm_hosp_30dmort)
train$hosp_30dmort <- as.numeric(as.character(train$hosp_30dmort)) 
test$hosp_30dmort <- as.numeric(as.character(test$hosp_30dmort)) 
SVM_pred <- as.numeric(a$svm_pred)


#XGBoost model
train$hosp_30dmort <- as.numeric(as.character(train$hosp_30dmort))#确保应变量为数值类型
test$hosp_30dmort <- as.numeric(as.character(test$hosp_30dmort))#确保应变量为数值类型
library(xgboost)
train_matrix <- xgb.DMatrix(data.matrix(train[,-1]),label = train$hosp_30dmort)
test_matrix <- xgb.DMatrix(data.matrix(test[,-1]),label = test$hosp_30dmort)
params <- list(objective = "binary:logistic",eval_metric="logloss",max_depth=3,eta = 0.1,
               gamma=0.5,colsample_bytree = 1,min_child_weight=1,subsample=0.5)
watchlist <- list(train = train_matrix,val=test_matrix)
xgb_model <- xgb.train(params = params,data = train_matrix,nrounds = 125,watchlist = watchlist,
                       early_stopping_rounds = 10,print_every_n = 10,maximize = FALSE)
xgb_pred_prob <- predict(xgb_model,test_matrix)
xgb_pred <- ifelse(xgb_pred_prob > 0.5,1,0)
xgb_pred_factor <- factor(xgb_pred,levels = c(0,1))
test_hosp_30dmort_factor <- factor(test$hosp_30dmort,levels = c(0,1))
confusionMatrix(data = xgb_pred_factor,reference=test_hosp_30dmort_factor)


#RF model（随机森林）
library(randomForest)
set.seed(123)
train$hosp_30dmort <- as.factor(train$hosp_30dmort)
test$hosp_30dmort <- as.factor(test$hosp_30dmort)
rf_model <- randomForest(hosp_30dmort~.,data = train,ntree = 500,mtry = 6)
rf_pred <- predict(rf_model,newdata = test)
confusionMatrix(data = rf_pred,reference = test$hosp_30dmort)


#DT model(决策树)
library(rpart)
dt_model <- rpart(hosp_30dmort~.,data = train,method = "class")
dt_pred_prob <- predict(dt_model,newdata = test,type="prob")[,2]
dt_pred <- ifelse(dt_pred_prob > 0.5,1,0)
confusionMatrix(factor(dt_pred,levels = c("0","1")),test$hosp_30dmort)


#NB model(朴素贝叶斯)
library(e1071)
nb_model <- naiveBayes(hosp_30dmort~.,data = train)
nb_pred_prob <- predict(nb_model,newdata = test,type="raw")[,2]
nb_pred <- ifelse(nb_pred_prob > 0.5,1,0)
confusionMatrix(factor(nb_pred,levels = c("0","1")),test$hosp_30dmort)


#KNN model(K最近邻回归)
library(kknn)
knn_model <- kknn(hosp_30dmort~.,train,test,k=10,distance = 2,kernel = "rectangular")
knn_pred_prob <- predict(knn_model,newdata=test,type="prob")
knn_pred_prob <- knn_pred_prob[,"1"]
knn_pred_prob <- as.numeric(knn_pred_prob)
threshold <- 0.5
knn_pred <- ifelse(knn_pred_prob > threshold,1,0)
confusionMatrix(factor(knn_pred,levels = c("0","1")),test$hosp_30dmort)


#内部验证ROC图绘制(pROC)
library(pROC)
library(ggplot2)
test_hosp_30dmort <- test$hosp_30dmort
ML_ROC <- data.frame(test_hosp_30dmort,
                     lm_pred,LR_pred,
                     svm_pred,svm_pred_prob,
                     xgb_pred_prob,xgb_pred,
                     rf_pred,
                     dt_pred_prob,dt_pred,
                     nb_pred_prob,nb_pred,
                     knn_pred,knn_pred_prob
)
ML_ROC$test_hosp_30dmort <- as.numeric(as.character(ML_ROC$test_hosp_30dmort))
ML_ROC$svm_pred <- as.numeric(as.character(ML_ROC$svm_pred))
ML_ROC$rf_pred <- as.numeric(as.character(ML_ROC$rf_pred))
str(ML_ROC)
#指定结局变量并命名
LR_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$lm_pred);LR_ROC
SVM_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$svm_pred_prob);SVM_ROC
XGB_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$xgb_pred_prob);XGB_ROC
RF_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$rf_pred);RF_ROC
DT_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$dt_pred_prob);DT_ROC
NB_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$nb_pred_prob);NB_ROC
KNN_ROC <- roc(ML_ROC$test_hosp_30dmort,ML_ROC$knn_pred_prob);KNN_ROC

#输出AUC的95%置信区间
auc(LR_ROC)
auc(SVM_ROC)
auc(XGB_ROC)
auc(RF_ROC)
auc(DT_ROC)
auc(NB_ROC)
auc(KNN_ROC)
ci.auc(LR_ROC)
ci.auc(SVM_ROC)
ci.auc(XGB_ROC)
ci.auc(RF_ROC)
ci.auc(DT_ROC)
ci.auc(NB_ROC)
ci.auc(KNN_ROC)

#指定图像保存路径和文件名
#png("ML ROC curves.png",width = 800,height = 800)
jpeg(filename = "Figure 4B 外部验证ROC.jpeg", width = 6, height = 6, units = "in", res = 600)
#画出第一条曲线
plot.roc(SVM_ROC,
         max.auc.polygon=F,#填充整个图像
         smooth=F,#绘制不平滑曲线“
         #main="Comparison of different machine learning models of ROC curves",#添加标题
         col = "red",#曲线颜色为红色
         legacy.axes=T,#使横轴从0到1，表示为1-特异度
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.8
)
#逐步添加其他曲线
#plot.roc(SVM_ROC,
#         add = T,
#         col = "orange",
#         smooth = F,
#         lwd=2
#         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.75
#)
plot.roc(XGB_ROC,
         add = T,
         col = "black",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.70
)
plot.roc(RF_ROC,
         add = T,
         col = "pink",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.65
)
plot.roc(DT_ROC,
         add = T,
         col = "green",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.60
)
plot.roc(NB_ROC,
         add = T,
         col = "blue",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.55
)
plot.roc(KNN_ROC,
         add = T,
         col = "orange",
         smooth = F,
         lwd=2
         #print.auc=TRUE,print.auc.x=0.2,print.auc.y=0.50
)
#增加图例
legend("bottomright",
       legend = c("SVM（AUC=0.6825）","XGBoost（AUC=0.697）","RF（AUC=0.5516）","DT（AUC=0.6328）","NB（AUC=0.6608）","KNN（AUC=0.6187）"),
       col = c("red","black","pink","green","blue","orange"),
       lwd = 1,
       cex=0.7
)
dev.off()

#将测试集重新切换为内部验证集，并且将内部验证集重新跑一边，然后进行校准曲线与决策曲线分析
test <- mimic_sad_ml[-index,]
#校准曲线绘制
#使用rms包绘制
library(rms)
ML_ROC$test_hosp_30dmort <- as.factor(as.character(ML_ROC$test_hosp_30dmort))
fit_SVM <- lrm(test_hosp_30dmort~svm_pred_prob,data = ML_ROC,x=TRUE,y=TRUE)
fit_XGB <- lrm(test_hosp_30dmort~xgb_pred_prob,data = ML_ROC,x=TRUE,y=TRUE)
fit_RF <- lrm(test_hosp_30dmort~rf_pred,data = ML_ROC,x=TRUE,y=TRUE)
fit_DT <- lrm(test_hosp_30dmort~dt_pred_prob,data = ML_ROC,x=TRUE,y=TRUE)
fit_NB <- lrm(test_hosp_30dmort~xgb_pred_prob,data = ML_ROC,x=TRUE,y=TRUE)
fit_KNN <- lrm(test_hosp_30dmort~knn_pred_prob,data = ML_ROC,x=TRUE,y=TRUE)
cal_SVM <- calibrate(fit_SVM,method="boot",B=1000)
cal_XGB <- calibrate(fit_XGB,method="boot",B=1000)
cal_RF <- calibrate(fit_RF,method="boot",B=1000)
cal_DT <- calibrate(fit_DT,method="boot",B=1000)
cal_NB <- calibrate(fit_NB,method="boot",B=1000)
cal_KNN <- calibrate(fit_KNN,method="boot",B=1000)
jpeg(filename = "Figure 4C 校准曲线.jpeg", width = 6, height = 6, units = "in", res = 600)
plot(1,type="n",#绘制一幅空白图
     xlim=c(0,1),#设置X轴范围为0-1
     ylim=c(0,1),#设置y轴范围为0-1
     xaxs="i",#设置图片原点
     yaxs="i",#设置图片原点
     xlab="Predicted Probability",#设置x轴名称
     ylab="Observed Probability",#设置y轴名称
     legend=FALSE,#不显示图例
     subtitle=FALSE,#不显示副标题
     #设置坐标轴刻度线及名称的字体相对大小
     cex=1.2,
     cex.axis=1.2,
     cex.lab=1.2)
abline(0,1,col="grey",lty=2,lwd=1)
lines(cal_SVM[,c("predy","calibrated.corrected")],lty=1,lwd=1.5,col="red")
lines(cal_XGB[,c("predy","calibrated.corrected")],lty=1,lwd=1.5,col="black")
lines(cal_RF[,c("predy","calibrated.corrected")],lty=1,lwd=1.5,col="pink")
lines(cal_DT[,c("predy","calibrated.corrected")],lty=1,lwd=1.5,col="green")
lines(cal_NB[,c("predy","calibrated.corrected")],lty=1,lwd=1.5,col="blue")
lines(cal_KNN[,c("predy","calibrated.corrected")],lty=1,lwd=1.5,col="orange")
legend(0.01,0.95,
       c("SVM",
         "XGBoost",
         "RF",
         "DT",
         "NB",
         "KNN"),
       lty=c(1,1,1,1,1,1),
       lwd=c(1.5,1.5,1.5,1.5,1.5,1.5),
       col = c("red","black","pink","green","blue","orange"),
       bty = "n",#是否显示图例边框
       cex=1.1)
dev.off()

#决策曲线绘制
#使用dcurves包
library(dcurves)
#整理数据类型
ML_DCA_cur <- ML_ROC
ML_DCA_cur$test_hosp_30dmort <- as.numeric(as.character(ML_DCA_cur$test_hosp_30dmort))
ML_DCA_cur$lm_pred <- as.numeric(as.character(ML_DCA_cur$lm_pred))
ML_DCA_cur$svm_pred_prob <- as.numeric(as.character(ML_DCA_cur$svm_pred_prob))
ML_DCA_cur$xgb_pred_prob <- as.numeric(as.character(ML_DCA_cur$xgb_pred_prob))
ML_DCA_cur$rf_pred <- as.numeric(as.character(ML_DCA_cur$rf_pred))
ML_DCA_cur$dt_pred_prob <- as.numeric(as.character(ML_DCA_cur$dt_pred_prob))
ML_DCA_cur$nb_pred_prob <- as.numeric(as.character(ML_DCA_cur$nb_pred_prob))
ML_DCA_cur$knn_pred_prob <- as.numeric(as.character(ML_DCA_cur$knn_pred_prob))

jpeg(filename = "Figure 4D 决策曲线.jpeg", width = 7, height = 6, units = "in", res = 600)
#绘制DCA曲线
dcurves::dca(formula = test_hosp_30dmort~svm_pred_prob+xgb_pred_prob+rf_pred+dt_pred_prob+nb_pred_prob+knn_pred_prob,
             label = list(svm_pred_prob="SVM",xgb_pred_prob="XGBoost",rf_pred="RF",dt_pred_prob="DT",nb_pred_prob="NB",knn_pred_prob="KNN"),
             data = ML_DCA_cur) %>% 
  plot(smooth=TRUE)+
  ggplot2::labs(x="Threshold Probability")+
  theme_bw()
dev.off()

##SHAP值进行模型解释
library(shapviz)
library(SHAPforxgboost)
library(tidyverse)
library(caret)
getSelectedAttributes(boruta_output)#获取确认的特征名称
train$hosp_30dmort <- as.numeric(as.character(train$hosp_30dmort))
train$age <- as.numeric(as.character(train$age))
train$weight <- as.numeric(as.character(train$weight))
train$temperature <- as.numeric(as.character(train$temperature))
train$heart_rate <- as.numeric(as.character(train$heart_rate))
train$resp_rate <- as.numeric(as.character(train$resp_rate))
train$spo2 <- as.numeric(as.character(train$spo2))
train$mbp <- as.numeric(as.character(train$mbp))
train$wbc <- as.numeric(as.character(train$wbc))
train$hemoglobin <- as.numeric(as.character(train$hemoglobin))
train$bun <- as.numeric(as.character(train$bun))
train$Cl <- as.numeric(as.character(train$Cl))
train$K <- as.numeric(as.character(train$K))
train$P <- as.numeric(as.character(train$P))
train$inr <- as.numeric(as.character(train$inr))
train$ptt <- as.numeric(as.character(train$ptt))
train$bicarbonate <- as.numeric(as.character(train$bicarbonate))
train$aniongap <- as.numeric(as.character(train$aniongap))
train$crrt <- as.numeric(as.character(train$crrt))
train$vaso <- as.numeric(as.character(train$vaso))
train$seda <- as.numeric(as.character(train$seda))
train$sofa_score <- as.numeric(as.character(train$sofa_score))
train$ckd <- as.numeric(as.character(train$ckd))
train$aki <- as.numeric(as.character(train$aki))
train$stroke <- as.numeric(as.character(train$stroke))
#train中括号里面的数据根据变量数量调整，24个变量加应变量，就是从第2个到25个
shap_values <- shapviz::shapviz(xgb_model,X_pred = as.matrix(train[2:25]))
jpeg(filename = "Figure 5A SHAP柱状图.jpeg", width = 6, height = 6, units = "in", res = 600)
sv_importance(shap_values,max_display = Inf,kind = "bar",show_numbers = T)+theme_bw()
dev.off()
jpeg(filename = "Figure 5B SHAP蜂群图.jpeg", width = 6, height = 6,units = "in", res = 600)
sv_importance(shap_values,max_display = Inf,kind = "beeswarm",show_numbers = T)+theme_bw()
dev.off()
                                                    
#绘制部分依赖图
#library(xgboost)
#library(pdp)
#library(ggplot2)
#library(shapviz)
#计算SHAP值
#shap_data <- shapviz(xgb_model,X_pred = as.matrix(train[2:25]))
#绘制部分依赖图
#jpeg(filename = "Figure 6A.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "bun",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6B.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "age",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6C.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "pt",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6D.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "platelet",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6E.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "resp_rate",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6F.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "ptt",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6G.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "temperature",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6H.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "aki",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(method="lm",se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6I.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "heart_rate",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6J.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "sbp",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6K.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "hemoglobin",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6L.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "Cl",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6M.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "weight",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6N.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "glu",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(se = FALSE,color="red")
#dev.off()

#jpeg(filename = "Figure 6O.jpeg", width = 6, height = 6,units = "in", res = 600)
#sv_dependence(shap_data,v = "vaso",alpha = 0.5,size = 1.5,color_var = NULL,color = "#3b528b",jitter_width = 0.05)+ 
#  geom_smooth(method="lm",se = FALSE,color="red")
#dev.off()


#将机器学习模型效能表格(热图)
# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

#内部验证集
# 创建数据框
data_internal_test <- data.frame(
  Model = c("SVM", "XGBoost", "RF", "DT", "NB", "KNN"),
  AUROC = c(0.6996,0.7308,0.5295,0.6388,0.6873,0.6645),
  Sensitivity = c(0.9804,0.9671,0.9859,0.9733,0.8447,0.9812),
  Specificity = c(0.1082,0.1228,0.0731,0.0585,0.3684,0.0614),
  PPV = c(0.8039,0.8043,0.7986,0.7940,0.8329,0.7958),
  NPV = c(0.5968,0.5000,0.5814,0.3704,0.3889,0.4667),
  Accuracy = c(0.7959,0.7885,0.7928,0.7798,0.7440,0.7866),
  Kappa = c(0.1265,0.1241,0.0867,0.0460,0.2173,0.0624)
)

# 转换为长格式
data_long_internal_test <- data_internal_test %>%
  pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")

# 绘制热图
heatmap_internal_test <- ggplot(data_long_internal_test, aes(x = Metric, y = Model, fill = Value)) +
                         geom_tile(color = "white") +
                         geom_text(aes(label = sprintf("%.4f", Value)), color = "black", size = 3) +
                         scale_fill_distiller(palette = "RdYlBu") + # 红-黄-蓝渐变色 
                         labs(
                         title = "Model performance on the internal validation set",
                         x = "Model performance metrics",
                         y = "Model",
                         fill = "score"
                         ) +
                         theme_minimal() +
                         theme(
                         axis.text.x = element_text(angle = 45, hjust = 1),
                         plot.title = element_text(hjust = 0.5, face = "bold")
                         )+
                         coord_flip() # 旋转坐标轴，使模型名称在底部，性能指标在左侧
jpeg(filename = "内部验证集模型性能热图.jpeg", width = 7, height = 6,units = "in", res = 600)
print(heatmap_internal_test)
dev.off()

#外部验证集
# 创建数据框
data_external_test <- data.frame(
  Model = c("SVM", "XGBoost", "RF", "DT", "NB", "KNN"),
  AUROC = c(0.6825,0.6970,0.5516,0.6328,0.6608,0.6187),
  Sensitivity = c(0.9500,0.9690,0.9714,0.9571,0.5929,0.9691),
  Specificity = c(0.1538,0.1758,0.1319,0.0550,0.6044,0.0879),
  PPV = c(0.8382,0.8444,0.8378,0.8238,0.8737,0.8306),
  NPV = c(0.4000,0.5517,0.5000,0.2174,0.2434,0.3810),
  Accuracy = c(0.8082,0.8278,0.8219,0.7965,0.5949,0.8121),
  Kappa = c(0.1368,0.1976,0.1452,0.0171,0.1248,0.0815)
)

# 转换为长格式
data_long_external_test <- data_external_test %>%
  pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")

# 绘制热图
heatmap_external_test <- ggplot(data_long_external_test, aes(x = Metric, y = Model, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.4f", Value)), color = "black", size = 3) +
  scale_fill_distiller(palette = "RdYlBu") + # 红-黄-蓝渐变色 
  labs(
    title = "Model performance on the external validation set",
    x = "Model performance metrics",
    y = "Model",
    fill = "score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+
  coord_flip() # 旋转坐标轴，使模型名称在底部，性能指标在左侧
jpeg(filename = "外部验证集模型性能热图.jpeg", width = 7, height = 6,units = "in", res = 600)
print(heatmap_external_test)
dev.off()

#对MIMIC-IV数据进行COX回归分析
library(tidyverse)#数据处理包
library(haven)#加载dta格式数据
mimic_data_cox <- read_dta("MIMIC-IV.dta")

#筛选SAD患者
mimic_sad_cox <- mimic_data_cox %>% filter(sad==1)

#观察缺失值
library(naniar)
missing_summary_mimic <- miss_var_summary(mimic_data_cox)#显示各列缺失值占比
print(missing_summary_mimic)
gg_miss_var(mimic_data_cox) + theme_bw()#缺失值可视化

#生成患者住院 30天院内死亡率生存时间数据
library(survival)
library(survminer)
mimic_km_cox <- mimic_data_cox %>% select(hosp_mort,hospstay,everything())
mimic_km_cox <- mimic_km_cox %>% mutate(hospstay_d = round(hospstay/60))
mimic_km_cox <- mimic_km_cox %>% mutate(hospstay_30d = ifelse(hospstay_d >= 30,30,hospstay_d))
mimic_km_cox <- mimic_km_cox %>% mutate(hosp_30dmort = case_when(
  hospstay_d > 30 & hosp_mort == 1 ~ 0,
  hospstay_d > 30 & hosp_mort == 0 ~ 0,
  hospstay_d <= 30 & hosp_mort == 1 ~ 1,
  hospstay_d <= 30 & hosp_mort == 0 ~ 0
))#建立新的30天院内死亡率数据

#整理生存分析数据
mimic_km_cox <- mimic_km_cox %>% select(-c(race,first_careunit,deliriumtime,sepsistime,icu28dmort,hosp_mort,hospstay,hospstay_d,icustay,stay_id)) %>% select(hospstay_30d,hosp_30dmort,everything())

#使用MICE包对生存分析数据进行多重插补
library(mice)
set.seed(123)
mimic_km_cox <- mice(mimic_km_cox,m=5,maxit = 10,method = "pmm",seed = 123,print = TRUE)#设置种子，输出过程
mimic_km_cox <- mice::complete(mimic_km_cox)
mimic_km_cox <- mimic_km_cox %>% filter(sad==1)
mimic_km_cox <- mimic_km_cox[,-41]
#COX分析（KIMI版）
library(survival)
library(broom)
#单因素COX回归分析
#初始化一个向量来存储显著变量的名称
significant_vars <- c()
results_mimic_cox_list <- list()
results_cox_ci_list <- list()
#对第3到42列的每个变量进行单因素 COX 回归分析
for (i in 3:42) {
  # 提取变量名
  var_name <- colnames(mimic_km_cox)[i]
  # 拟合单因素 COX 回归模型
  cox_model <- coxph(Surv(mimic_km_cox[,1], mimic_km_cox[,2]) ~ mimic_km_cox[,i], data = mimic_km_cox)
  # 存储结果
  results_mimic_cox_list[[colnames(mimic_km_cox)[i]]] <- summary(cox_model)
  # 提取回归系数和置信区间
  coefficients <- summary(cox_model)$coefficients
  ci <- exp(confint(cox_model))  # 置信区间指数运算
  
  # 存储结果
  results_cox_ci_list[[colnames(mimic_km_cox)[i]]] <- data.frame(
    OR = exp(coefficients[, 1]),  # OR 值
    Lower_CI = ci[, 1],           # 置信区间下限
    Upper_CI = ci[, 2],           # 置信区间上限
    P_value = coefficients[, "Pr(>|z|)"]   # p 值
  )
  # 提取 P 值
  p_value <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
  # 判断 P 值是否小于 0.05
  if (p_value < 0.05) {
    significant_vars <- c(significant_vars, var_name)
  }
  
}

# 查看结果
for (var_name in names(results_cox_ci_list)) {
  cat("变量:", var_name, "\n")
  print(results_cox_ci_list[[var_name]])
  cat("\n")
}

#多因素COX分析
# 创建多因素 COX 回归的公式
multivariate_formula <- as.formula(paste("Surv(mimic_km_cox[,1], mimic_km_cox[,2]) ~", paste(significant_vars, collapse = " + ")))
# 拟合多因素 COX 回归模型
multivariate_cox_model <- coxph(multivariate_formula, data = mimic_km_cox)
# 查看多因素 COX 回归模型的结果
summary(multivariate_cox_model)

#提取模型结果
# 提取模型的系数、P值等信息
model_results <- summary(multivariate_cox_model)$coefficients
# 将模型结果转换为数据框
model_df <- as.data.frame(model_results)
# 获取模型中的变量名
var_names <- rownames(model_df)
model_df$Variable <- var_names
# 选择要导出的列
export_df <- model_df[, c("coef", "exp(coef)","Pr(>|z|)", "Variable")]
# 重命名列
colnames(export_df) <- c("Coefficient", "Hazard Ratio", "P Value", "Variable")
# 导出为CSV文件
write.csv(export_df, "多重COX回归 multivariate_cox_results.csv", row.names = FALSE)



#对SOFA、SAPS II、APS III评分预测SAD患者发生30天死亡进行ROC曲线分析
# 安装并加载必要的包
library(pROC)
library(tidyverse)
#给MIMIC-IV数据库纳入研究患者加入SAPSII与APSIII评分
library(tidyverse)
#导入SAPSII与SAPSIII评分数据
SAPSII <- read_csv("score_sapsii.csv")
APSIII <- read_csv("score_apsiii.csv")
mimic_sad <- mimic_data %>% filter(sad==1)
#建立MIMIC-IV新的30天院内死亡率数据
mimic_sad <- mimic_sad %>% mutate(hospstay_d = round(hospstay/60))
mimic_sad <- mimic_sad %>% mutate(hospstay_30d = ifelse(hospstay_d >= 30,30,hospstay_d))
mimic_sad <- mimic_sad %>% mutate(hosp_30dmort = case_when(
  hospstay_d > 30 & hosp_mort == 1 ~ 0,
  hospstay_d > 30 & hosp_mort == 0 ~ 0,
  hospstay_d <= 30 & hosp_mort == 1 ~ 1,
  hospstay_d <= 30 & hosp_mort == 0 ~ 0
))
#建立新的数据框内容包括SOFA评分、SAPSII评分、APSIII评分、30天死亡
sofa_sad <- mimic_sad %>% select(stay_id,hosp_30dmort,sofa_score) 
sapsii_sofa_sad <- left_join(sofa_sad,SAPSII,by="stay_id")
apsiii_sapsii_sofa_sad <- left_join(sapsii_sofa_sad,APSIII,by="stay_id")
multiscore_sad <- apsiii_sapsii_sofa_sad %>% select(hosp_30dmort,sofa_score,sapsii,apsiii)
#判断数据框中是否存在缺失值
any(is.na(multiscore_sad))
#存在缺失值则计算每一列的缺失值比例
sapply(multiscore_sad, function(x) sum(is.na(x)) / nrow(multiscore_sad))
# 查找缺失值的位置
print(which(is.na(multiscore_sad), arr.ind = TRUE))
#对缺失值进行多重插补
library(mice)
multiscore_sad <- mice(multiscore_sad,m=5,maxit = 10,method = "pmm",seed = 123,print = TRUE)#设置种子，不输出过程
multiscore_sad <- mice::complete(multiscore_sad)

# 确保30天死亡是因子变量，其他各类评分是数值变量
multiscore_sad$hosp_30dmort <- as.factor(multiscore_sad$hosp_30dmort)
multiscore_sad$sofa_score <- as.numeric(multiscore_sad$sofa_score)
multiscore_sad$sapsii <- as.numeric(multiscore_sad$sapsii)
multiscore_sad$apsiii <- as.numeric(multiscore_sad$apsiii)

library(caret)
set.seed(666)
index <- createDataPartition(multiscore_sad$hosp_30dmort,p=0.7,list = F)
train_multi <- multiscore_sad[index,]
test_multi <- multiscore_sad[-index,]

# 构建单个评分逻辑回归模型
#SOFA评分
model_LR_sofa <- glm(hosp_30dmort ~ sofa_score, data = train_multi, family = binomial(link = "logit"))
summary(model_LR_sofa)
LR_sofa_pred <- predict(model_LR_sofa,test_multi,type="response")
threshold <- 0.5
predictions_binary_sofa <- ifelse(LR_sofa_pred>threshold,1,0)
confusionMatrix(as.factor(predictions_binary_sofa),as.factor(test_multi$hosp_30dmort))

#SAPS II评分
model_LR_sapsii <- glm(hosp_30dmort ~ sapsii,data = train_multi,family = binomial(link = "logit"))
summary(model_LR_sapsii)
LR_sapsii_pred <- predict(model_LR_sapsii,test_multi,type="response")
threshold <- 0.5
predictions_binary_sapsii <- ifelse(LR_sapsii_pred>threshold,1,0)
confusionMatrix(as.factor(predictions_binary_sapsii),as.factor(test_multi$hosp_30dmort))

#APS III评分
model_LR_apsiii <- glm(hosp_30dmort~apsiii,data = train_multi,family = binomial(link = "logit"))
summary(model_LR_apsiii)
LR_apsiii_pred <- predict(model_LR_apsiii,test_multi,type="response")
threshold <- 0.5
predictions_binary_apsiii <- ifelse(LR_apsiii_pred>threshold,1,0)
confusionMatrix(as.factor(predictions_binary_apsiii),as.factor(test_multi$hosp_30dmort))

#构建传统联合评分逻辑回归模型
model_LR_multi <- glm(hosp_30dmort~sofa_score + sapsii + apsiii,data = train_multi,family = binomial(link = "logit"))
summary(model_LR_multi)
LR_multi_pred <- predict(model_LR_multi,test_multi,type="response")
threshold <- 0.5
predictions_binary_multi <- ifelse(LR_multi_pred>threshold,1,0)
confusionMatrix(as.factor(predictions_binary_multi),as.factor(test_multi$hosp_30dmort))


# 计算ROC曲线和AUC
#SOFA评分
roc_curve_sofa <- roc(test_multi$hosp_30dmort,LR_sofa_pred)
auc(roc_curve_sofa)
ci.auc(roc_curve_sofa)
#SAPS II评分
roc_curve_sapsii <- roc(test_multi$hosp_30dmort,LR_sapsii_pred)
auc(roc_curve_sapsii)
ci.auc(roc_curve_sapsii)
#APS III评分
roc_curve_apsiii <- roc(test_multi$hosp_30dmort,LR_apsiii_pred)
auc(roc_curve_apsiii)
ci.auc(roc_curve_apsiii)
#Multi score评分
roc_curve_multi <- roc(test_multi$hosp_30dmort,LR_multi_pred)
auc(roc_curve_multi)
ci.auc(roc_curve_multi)

jpeg(filename = "Figure 7.jpeg", width = 2500, height = 2500, res = 400)
#画出第一条曲线
plot.roc(roc_curve_sofa,
         max.auc.polygon=F,#填充整个图像
         smooth=F,#绘制不平滑曲线“
         #main="Comparison of different machine learning models of ROC curves",#添加标题
         col = "darkgreen",#曲线颜色为绿色
         legacy.axes=T,#使横轴从0到1，表示为1-特异度
         lwd=2,
         print.auc=F,print.auc.x=0.2,print.auc.y=0.8
)
#逐步添加其他曲线
plot.roc(roc_curve_sapsii,
         add = T,
         col = "orange",
         smooth = F,
         lwd=2,
         print.auc=F,print.auc.x=0.2,print.auc.y=0.75)
plot.roc(roc_curve_apsiii,
         add = T,
         col = "black",
         smooth = F,
         lwd=2,
         print.auc=F,print.auc.x=0.2,print.auc.y=0.70)
plot.roc(roc_curve_multi,
         add = T,
         col = "blue",
         smooth = F,
         lwd=2,
         print.auc=F,print.auc.x=0.2,print.auc.y=0.65)
plot.roc(XGB_ROC,
         add = T,
         col = "red",
         smooth = F,
         lwd=2,
         print.auc=F,print.auc.x=0.2,print.auc.y=0.60)

#增加图例
legend("bottomright",
       legend = c("SOFA (AUC=0.586)","SAPS II (AUC=0.700)","APS III (AUC=0.696)","Combined score(AUC=0.715)","XGBoost (AUC=0.731)"),
       col = c("darkgreen","orange","black","blue","red"),
       lwd = 1,
       cex=0.9
)
dev.off()
