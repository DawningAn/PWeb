library(knitr)
library(ape)
library(MASS)
library(fields)
library(Hmsc)
library(corrplot)

# Set seed for reproducibility
# 设置随机种子
set.seed(1)

# Number of species
ns = 50

# Generate a coalescent phylogeny
# 生成系统发育树
phy = ape::rcoal(n = ns, tip.label = sprintf('species_%.3d', 1:ns), br = "coalescent")
# Plot the phylogeny
# 绘制系统发育树
plot(phy, show.tip.label = FALSE, no.margin = TRUE)

# Compute the variance-covariance matrix# based on the phylogeny# 计算协方差矩阵
C = vcv(phy, model = "Brownian", corr = TRUE)
spnames = colnames(C)

# Simulate traits
# 模拟物种性状
traits = matrix(NA, ncol = 2, nrow = ns)
for (i in 1:2) {
  traits[, i] = mvrnorm(n = 1, mu = rep(0, ns), Sigma = C)
}
rownames(traits) = spnames
colnames(traits) = c("habitat.use", "thermal.optimum")
traits = as.data.frame(traits)

# Set graphics parameters
# 设置图形参数并绘制图像
par(fig = c(0, 0.6, 0, 0.8), mar = c(6, 0, 2, 0))
plot(phy, show.tip.label = FALSE)
par(fig = c(0.6, 0.9, 0.025, 0.775), mar = c(6, 0, 2, 0), new = TRUE)

plot.new()
image.plot(t(traits), axes = FALSE, legend.width = 3, legend.shrink = 1, col = colorRampPalette(c("blue", "white", "red"))(200))
# Simulation of environmental data
# 模拟环境变量
n = 200
habitat = factor(sample(x = c("forest", "open"), size = n, replace = TRUE))
climate = rnorm(n)
nc = 4
mu = matrix(0, nrow = nc, ncol = ns)

# Expected niche calculations
# 预期生态位计算
mu[1,] = -traits$thermal.optimum^2 / 4
-traits$habitat.use
mu[2,] = 2 * traits$habitat.use
mu[3,] = traits$thermal.optimum / 2
mu[4,] = -1 / 4
beta = mu + 0.25 * matrix(rnorm(n = ns * nc), ncol = ns)
X = cbind(rep(1, ns), as.numeric(habitat == "forest"), climate, climate * climate)
L = X %*% beta
Y = L + mvrnorm(n = n, mu = rep(0, ns), Sigma = diag(ns))
colnames(Y) = spnames

# Data preparation for HMSC
# HMSC数据准备
XData = data.frame(climate = climate, habitat = habitat)
XFormula = ~habitat + poly(climate, degree = 2, raw = TRUE)
TrFormula = ~habitat.use + thermal.optimum
studyDesign = data.frame(sample = sprintf('sample_%.3d', 1:n), stringsAsFactors = TRUE)
rL = HmscRandomLevel(units = studyDesign$sample)
rL$nfMax = 15
rL
# Model specification
m = Hmsc(Y = Y, XData = XData, XFormula = XFormula,
         TrData = traits, TrFormula = TrFormula,
         phyloTree = phy, studyDesign = studyDesign, ranLevels = list(sample = rL))
Y
# MCMC settings
nChains = 2
test.run = FALSE
if (test.run) {
  thin = 1
  samples = 100
  transient = 50
} else {
  thin = 10
  samples = 1000
  transient = 500
}
verbose = 0

# Run MCMC
m = sampleMcmc(m, thin = thin, samples = samples,
               transient = transient, nChains = nChains, nParallel = 1, verbose = verbose) # Set nParallel=1 to avoid port issues

######## 检验马尔可夫链蒙特卡罗（MCMC）模型收敛性
# 将模型结果转换为Coda对象，以便进行统计分析
mpost = convertToCodaObject(m)

# 设置图形参数，准备绘制六个直方图
par(mfrow = c(3, 2))

# 计算Beta参数的有效样本大小，并绘制其直方图
ess.beta = effectiveSize(mpost$Beta)
hist(ess.beta)

# 计算Beta参数的Gelman-Rubin诊断，并绘制其直方图
psrf.beta = gelman.diag(mpost$Beta,
                        multivariate = FALSE)$psrf
hist(psrf.beta)

# 计算Gamma参数的有效样本大小，并绘制其直方图
ess.gamma = effectiveSize(mpost$Gamma)
hist(ess.gamma)

# 计算Gamma参数的Gelman-Rubin诊断，并绘制其直方图
psrf.gamma = gelman.diag(mpost$Gamma,
                         multivariate = FALSE)$psrf
hist(psrf.gamma)

# 从所有可能的物种对中随机抽样100个，用于进一步分析
sppairs = matrix(sample(x = 1:ns^2, size = 100))

# 提取并过滤Omega参数的后验样本
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)) {
  tmp[[chain]] = tmp[[chain]][, sppairs]
}

# 计算过滤后的Omega参数的有效样本大小，并绘制其直方图
ess.omega = effectiveSize(tmp)
hist(ess.omega)

# 计算过滤后的Omega参数的Gelman-Rubin诊断，并绘制其直方图
psrf.omega = gelman.diag(tmp,
                         multivariate = FALSE)$psrf
hist(psrf.omega)

# 打印Rho参数的有效样本大小
print("ess.rho:")
effectiveSize(mpost$Rho)

# 打印Rho参数的Gelman-Rubin诊断
print("psrf.rho:")
gelman.diag(mpost$Rho)$psrf


######## 评估模型拟合度和进行方差分解
# 计算模型的预测值
preds = computePredictedValues(m)

# 评估模型拟合度，包括计算R平方值，并绘制其直方图
MF = evaluateModelFit(hM = m, predY = preds)
hist(MF$R2, xlim = c(0, 1),
     main = paste0("Mean = ", round(mean(MF$R2), 2)))
# 显示模型中使用的环境变量
head(m$X)

# 计算方差分解，分析环境因子（如栖息地和气候）# 对模型解释的贡献
VP = computeVariancePartitioning(m,
                                 group = c(1, 1, 2, 2),
                                 groupnames = c("habitat", "climate"))
# 绘制方差分解的结果
plotVariancePartitioning(m, VP = VP)

# 使用kable函数以表格形式展示Beta参数的方差分解结果
kable(VP$R2T$Beta)

# 输出Y变量的方差分解结果
VP$R2T$Y

# 获取Beta参数的后验估计，并绘制其支持水平图
postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support",
         plotTree = TRUE, supportLevel = 0.95, split = .4, spNamesNumbers = c(F, F))
# 获取Gamma参数的后验估计，并绘制其支持水平图
postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post = postGamma, param = "Support", supportLevel = 0.95)
# 计算物种间的关联性，并根据支持水平进行筛选，绘制关联图
OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support > supportLevel)
  + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(200), tl.cex = .6, tl.col = "black",
         title = paste("random effect level:", m$rLNames[1]), mar = c(0, 0, 1, 0))
# 输出Rho参数的统计摘要
summary(mpost$Rho)


######## 构建和分析模型梯度，以探索特定环境变量## （如气候和栖息地）对物种分布的影响
# 构建以气候为焦点变量的梯度，其中栖息地被视为非焦点变量
Gradient = constructGradient(m, focalVariable = "climate", non.focalVariables = list("habitat" = list(3, "open")))
# 查看新构建的环境数据
Gradient$XDataNew

# 根据新的梯度数据进行预测，并计算预期值
predY = predict(m, XData = Gradient$XDataNew, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = TRUE)

# 绘制关于多样性指标S的梯度图，展示数据点
plotGradient(m, Gradient, pred = predY, measure = "S", showData = TRUE)
# 绘制关于生态响应Y的梯度图，指定显示第一个指标，展示数据点
plotGradient(m, Gradient, pred = predY, measure = "Y", index = 1, showData = TRUE)
# 绘制关于性状T的梯度图，指定第三个指标，展示数据点
plotGradient(m, Gradient, pred = predY, measure = "T", index = 3, showData = TRUE)
# 构建以栖息地为焦点变量的梯度，气候作为非焦点变量
Gradient = constructGradient(m, focalVariable = "habitat", non.focalVariables = list("climate" = list(1)))

# 查看新构建的环境数据
Gradient$XDataNew

# 根据新的梯度数据进行预测
predY = predict(m, XData = Gradient$XDataNew, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = TRUE)

# 绘制关于生态响应Y的梯度图，# 选择栖息地使用最高的指标，展示数据点，# 轻微调整数据点位置以避免重叠plotGradient(m, Gradient, pred=predY,   measure="Y", index=which.max(m$TrData$habitat.use),  showData = TRUE, jigger = 0.2)

# 绘制关于性状T的梯度图，选择第二个指标，# 展示数据点，轻微调整数据点位置
plotGradient(m, Gradient, pred = predY,
             measure = "T", index = 2, showData = TRUE, jigger = 0.2)
######## 错误指定模型的HMSC分析
# HMSC analyses of misspecified models
## Missing environmental covariate
## 缺少环境协变量
# 定义环境变量的公式，这里使用了气候数据的二次多项式形式
XFormula.1 = ~poly(climate, degree = 2, raw = TRUE)

# 创建HMSC模型对象，指定Y（响应变量数据），# XData（环境变量数据），XFormula（环境变量公式），# TrData（性状数据），TrFormula（性状公式），# phyloTree（系统发育树），studyDesign（研究设计），# ranLevels（随机效应层级）
ma50 = Hmsc(Y = Y, XData = XData, XFormula = XFormula.1,
            TrData = traits, TrFormula = TrFormula,
            phyloTree = phy,
            studyDesign = studyDesign, ranLevels = list(sample = rL))
# 对HMSC模型进行MCMC采样，指定采样间隔（thin），# 样本数（samples），过渡期（transient），# 链的数量（nChains），并行数量（nParallel），# 是否显示详细信息（verbose）
ma50 = sampleMcmc(ma50, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = 1, verbose = verbose)
# 计算方差分解，评估'climate'变量在模型中的贡献
VP = computeVariancePartitioning(ma50, group = c(1, 1, 1), groupnames = c("climate"))
# 绘制方差分解结果
plotVariancePartitioning(ma50, VP = VP)

# 计算物种间的关联性，并基于支持水平进行过滤
OmegaCor = computeAssociations(ma50)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support > supportLevel)
  + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean
# 使用corrplot绘制物种间关联性的热图，# 使用从蓝到红的颜色渐变
corrplot(toPlot, method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = .6, tl.col = "black",
         title = paste("random effect level:", ma50$rLNames[1]), mar = c(0, 0, 1, 0))
######## 在缺少某些性状数据时对模型分析的影响
## Missing traits

# 定义性状数据的公式，这里使用了栖息地使用情况作为性状变量
TrFormula.1 = ~habitat.use

# 创建HMSC模型对象，指定Y（响应变量数据），# XData（环境变量数据），XFormula（环境变量公式），# TrData（性状数据），TrFormula（性状数据公式），# phyloTree（系统发育树），studyDesign（研究设计），# ranLevels（随机效应层级）
m = Hmsc(Y = Y, XData = XData, XFormula = XFormula,
         TrData = traits, TrFormula = TrFormula.1,
         phyloTree = phy,
         studyDesign = studyDesign, ranLevels = list(sample = rL))
# 对HMSC模型进行MCMC采样，指定采样间隔（thin），# 样本数（samples），过渡期（transient），# 链的数量（nChains），并行数量（nParallel），# 是否显示详细信息（verbose）
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = 1, verbose = verbose)
# 计算方差分解，评估'habitat'和'climate'变量在模型中的贡献
VP = computeVariancePartitioning(m, group = c(1, 1, 2, 2), groupnames = c("habitat", "climate"))
# 绘制方差分解结果
plotVariancePartitioning(m, VP = VP)

# 使用kable函数以表格形式展示Beta参数的方差分解结果
kable(VP$R2T$Beta)

# 输出Y变量的方差分解结果
VP$R2T$Y

# 将模型结果转换为Coda对象，以便进行统计分析
mpost = convertToCodaObject(m)

# 输出Rho参数的统计摘要
summary(mpost$Rho)

########改变先验分布来影响物种负载的模型，# 并进行相关的统计分析和图形显示# Changing prior distribution for# the species loadings# 定义环境变量公式，这里使用了气候数据的二次多项式形式
XFormula.1 = ~poly(climate, degree = 2, raw = TRUE)

# 设置随机效应层级的先验分布参数
rL = setPriors(rL, a1 = 5, a2 = 5)

# 显示随机效应层级的结构
str(rL)

# 创建HMSC模型对象，并使用先前定义的先验分布参数
ma5 = Hmsc(Y = Y, XData = XData, XFormula = XFormula.1,
           TrData = traits, TrFormula = TrFormula,
           phyloTree = phy,
           studyDesign = studyDesign, ranLevels = list(sample = rL))
# 对HMSC模型进行MCMC采样
ma5 = sampleMcmc(ma5, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = 1, verbose = verbose)
# 更改先验分布参数
rL = HmscRandomLevel(units = studyDesign$sample)
rL = setPriors(rL, a1 = 500, a2 = 500)
str(rL)

# 创建新的HMSC模型对象，使用更新的先验分布参数
ma500 = Hmsc(Y = Y, XData = XData, XFormula = XFormula.1, TrData = traits, TrFormula = TrFormula, phyloTree = phy,
             studyDesign = studyDesign, ranLevels = list(sample = rL))
# 对更新的HMSC模型进行MCMC采样
ma500 = sampleMcmc(ma500, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = 1, verbose = verbose)
# 设置图形显示参数，准备显示三个关联图
par(mfrow = c(1, 3), mar = c(0, 0, 0, 0))

# 计算物种间关联性并绘制关联图，对于不同的先验设置分别显示
# a1 = a2 = 5
OmegaCor = computeAssociations(ma5)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support > supportLevel)
  + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$meancorrplot(toPlot, method = "color",
                                                                                   col = colorRampPalette(c("blue", "white", "red"))(200), tl.cex = .6, tl.col = "black",
                                                                                   title = "a1 = a2 = 5", mar = c(0, 0, 1, 0))

# a1 = a2 = 50
OmegaCor = computeAssociations(ma50)
toPlot = ((OmegaCor[[1]]$support > supportLevel)
  + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$meancorrplot(toPlot, method = "color",
                                                                                   col = colorRampPalette(c("blue", "white", "red"))(200), tl.cex = .6, tl.col = "black",
                                                                                   title = "a1 = a2 = 50", mar = c(0, 0, 1, 0))

# a1 = a2 = 500
OmegaCor = computeAssociations(ma500)
toPlot = ((OmegaCor[[1]]$support > supportLevel)
  + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$meancorrplot(toPlot, method = "color",
                                                                                   col = colorRampPalette(c("blue", "white", "red"))(200), tl.cex = .6, tl.col = "black",
                                                                                   title = "a1 = a2 = 500", mar = c(0, 0, 1, 0))