# 加载必要的库
suppressPackageStartupMessages(suppressWarnings(library(getopt)))
suppressPackageStartupMessages(suppressWarnings(library(magrittr)))
suppressPackageStartupMessages(suppressWarnings(library(ropls)))
suppressPackageStartupMessages(suppressWarnings(library(metabolomics)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(reshape2)))
suppressPackageStartupMessages(suppressWarnings(library(readr)))

# 定义命令行参数
commands = matrix(c(
  'help', 'h', 0, 'logical','Print this help file and exit.',
  'inputDir', 'i', 1, 'character','Directory of input data.',
  'outputDir', 'o', 1, 'character','Directory for output data.',
  'dataName', 'd', 1, 'character','Information of data file.',
  'sampleInfoName', 'a', 1, 'character','Name of Sample file.',
  'group', 'g', 1, 'character','The name of column that contains groups information in sampleInfo file.',
  'method', 'm', 1, 'character','PLS or OPLS.',
  'scaleC', 'S', 1, 'character','One of "none", "center", "pareto" and "standard".'
), byrow=TRUE, ncol=5)

opt = getopt(commands)

# 帮助信息
if (!is.null(opt$help)) {
  cat(paste(getopt(commands, usage = TRUE), "\n"))
  q(status=1)
}

# 参数检查
check_required_param <- function(param_name, param_value, error_message) {
  if (is.null(param_value)) {
    stop(paste0("Please set the parameter '", param_name, "': ", error_message))
  }
}

check_required_param("inputDir", opt$inputDir, "Input directory is required.")
check_required_param("outputDir", opt$outputDir, "Output directory is required.")
check_required_param("dataName", opt$dataName, "Data file name is required.")
check_required_param("sampleInfoName", opt$sampleInfoName, "Sample information file name is required.")
check_required_param("group", opt$group, "Group column name is required.")
check_required_param("method", opt$method, "Method (PLS or OPLS) is required.")
check_required_param("scaleC", opt$scaleC, "Scaling method is required.")

inputDir <- opt$inputDir
outputDir <- opt$outputDir
dataName <- opt$dataName
sampleInfoName <- opt$sampleInfoName
group <- opt$group
method <- opt$method
scaleC <- opt$scaleC

if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# 加载数据
dataFile <- paste0(inputDir, "/", dataName)
data <- read_csv(dataFile, col_types = cols(), progress = FALSE) %>% as.data.frame()
data <- data[!apply(data, 1, function(x) all(is.na(x))), ]
data <- data[, !apply(data, 2, function(x) all(is.na(x)))]
if (nrow(data) == 0) {
  stop("The data file is empty!")
}
rownames(data) <- data$compound.id

sampleInfoFile <- paste0(inputDir, "/", sampleInfoName)
sampleInfo <- read_csv(sampleInfoFile, col_types = cols()) %>% as.data.frame()
sampleInfo <- sampleInfo[!apply(sampleInfo, 1, function(x) all(is.na(x))), ]
sampleInfo <- sampleInfo[, !apply(sampleInfo, 2, function(x) all(is.na(x)))]
if (nrow(sampleInfo) == 0) {
  stop("The sample information file is empty!")
}
sampleInfo <- sampleInfo[sampleInfo$sample.names %in% colnames(data), ]
rownames(sampleInfo) <- sampleInfo$sample.names

# 转置数据
data <- data[, row.names(sampleInfo)] %>% t

# 分组信息
groups <- sampleInfo[, group]
groups <- as.factor(groups)
names(groups) <- row.names(sampleInfo)

# 执行 PLS-DA 或 OPLS-DA
setwd(outputDir)
if (method == "OPLS") {
  if (length(unique(groups)) != 2) {
    stop("OPLS-DA method is only suitable for 2-group data.")
  }
  oplsda <- opls(x = data, y = groups, scaleC = scaleC, fig.pdfC = "Rplots.pdf", info.txtC = paste0(outputDir, "/opls.txt"), predI = 1, orthoI = NA)
  system("mv Rplots.pdf opls.pdf 2>&1", intern = TRUE)
} else if (method == "PLS") {
  oplsda <- opls(x = data, y = groups, scaleC = scaleC, fig.pdfC = "Rplots.pdf", info.txtC = paste0(outputDir, "/pls.txt"), predI = NA, orthoI = 0)
  system("mv Rplots.pdf pls.pdf 2>&1", intern = TRUE)
} else {
  stop("Method must be 'PLS' or 'OPLS'.")
}

# 输出结果
write.csv(oplsda@scoreMN, paste0(outputDir, "/", method, ".score.csv"))
write.csv(oplsda@loadingMN, paste0(outputDir, "/", method, ".loading.csv"))

# 计算相关性
xModelMN <- oplsda@suppLs[["xModelMN"]]
yModelMN <- oplsda@suppLs[["yModelMN"]]
write.csv(xModelMN, paste0(outputDir, "/", method, ".xModelMN.csv"))
write.csv(yModelMN, paste0(outputDir, "/", method, ".yModelMN.csv"))

cxtCompMN <- cor(xModelMN, oplsda@scoreMN, use = "pairwise.complete.obs")
cytCompMN <- cor(yModelMN, oplsda@scoreMN, use = "pairwise.complete.obs")
write.csv(cxtCompMN, paste0(outputDir, "/", method, ".correlation.x_t.csv"))
write.csv(cytCompMN, paste0(outputDir, "/", method, ".correlation.y_t.csv"))

# 计算协方差
covx <- cov(xModelMN, oplsda@scoreMN, use = "pairwise.complete.obs")
covy <- cov(yModelMN, oplsda@scoreMN, use = "pairwise.complete.obs")
write.csv(covx, paste0(outputDir, "/", method, ".covariance.x_t.csv"))
write.csv(covy, paste0(outputDir, "/", method, ".covariance.y_t.csv"))

# S-plot 数据
splotdata <- cbind(pcorr1 = cxtCompMN[, 1], p1 = covx[, 1])
vip <- cbind(ortVIP = getVipVn(oplsda, orthoL = TRUE), VIP = getVipVn(oplsda))
rownames(vip) <- names(getVipVn(oplsda))
splotdata <- cbind(splotdata, vip)
write.csv(splotdata, paste0(outputDir, "/", method, ".splot_data.csv"))

# S-plot 图
vipCutoff <- 1
splotdata <- splotdata %>% mutate(color = ifelse(VIP > vipCutoff, TRUE, FALSE))
write.csv(splotdata, paste0(outputDir, "/", method, ".splot_data_plot.csv"))

pdf(paste0(outputDir, "/", method, ".splot.pdf"))
p <- ggplot(splotdata, aes(x = pcorr1, y = p1, color = color)) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 15)) +
  xlab("p") + ylab("pcorr") +
  geom_hline(aes(yintercept = 0), color = "grey") +
  geom_vline(aes(xintercept = 0), color = "grey") +
  geom_point() +
  scale_color_discrete(name = "", labels = c("VIP <= 1", "VIP > 1"))
print(p)
dev.off()

# 输出模型统计结果
if (method == "OPLS") {
  modBarDF <- oplsda@modelDF[!(rownames(oplsda@modelDF) %in% c("rot", "sum")), ]
} else {
  modBarDF <- oplsda@modelDF
}
write.csv(modBarDF, paste0(outputDir, "/", method, ".CV_results.csv"))

# 输出置换检验结果
permMN <- as.data.frame(oplsda@suppLs[["permMN"]])
write.csv(permMN, paste0(outputDir, "/", method, ".permutation_results.csv"), row.names = FALSE)

# 绘制置换检验图
permMN.plot <- permMN[, c(2, 3, 7)]
R2.y.mean <- mean(permMN.plot$`R2Y(cum)`)
Q2.y.mean <- mean(permMN.plot$`Q2(cum)`)
sim.x.mean <- mean(permMN.plot$sim)
permMN.plot.m <- melt(permMN.plot, id = 3)
title <- paste0("R2Y=", oplsda@summaryDF$pR2Y, ", ", "pQ2=", oplsda@summaryDF$pQ2)

pdf(paste0(outputDir, "/", method, ".permutation.pdf"))
p <- ggplot(permMN.plot.m, aes(x = sim, y = value, color = variable, shape = variable)) +
  geom_point() +
  geom_segment(aes(x = sim.x.mean, xend = permMN.plot[1, 3], y = R2.y.mean, yend = permMN.plot[1, 1]), color = "black") +
  geom_segment(aes(x = sim.x.mean, xend = permMN.plot[1, 3], y = Q2.y.mean, yend = permMN.plot[1, 2]), color = "black") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = 0) +
  labs(title = title, y = "", x = "Similarity ( y, y perm )") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(p)
dev.off()

# 删除临时文件
system("rm Rplots[1-5].pdf 2>&1", intern = TRUE)