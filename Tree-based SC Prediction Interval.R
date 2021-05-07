# Execute this file or source() it into an R session to launch

# =============================================================
# Tree-Based Synthetic Control Prediction Intervals


# =============================================================

# Import packages
if(!require(qwraps2)){install.packages('qwraps2');library(qwraps2)}
if(!require(pastecs)){install.packages('pastecs');library(pastecs)}
if(!require(Synth)){install.packages('Synth');library(Synth)}
if(!require(rsample)){install.packages('rsample');library(rsample)}
if(!require(randomForest)){install.packages('randomForest');library(randomForest)}
if(!require(rpart)){install.packages('rpart');library(rpart)}
if(!require(Rfast)){install.packages('Rfast');library(Rfast)}
if(!require(rpart.plot)){install.packages('rpart.plot');library(rpart.plot)}
if(!require(SCtools)){install.packages('SCtools');library(SCtools)}
if(!require(gsynth)){install.packages('gsynth');library(gsynth)}
if(!require(augsynth)){install.packages('augsynth');library(augsynth)}


# The dataset can be downloaded from the replication files of Cattaneo et al.(2021)
# https://github.com/mdcattaneo/replication-CFT_2021_JASA


# import data
setwd("~/desktop/replication-CFT_2021_JASA-master")

data <- read.csv("California.csv")
data$stateno <- as.numeric(as.factor(data$state))

# To test the result on Alabama, skip the following three lines.
data$stateno[data$state == "California"] <- 1
data$stateno[data$state == "Alabama"] <- 2
data$stateno[data$state == "Arkansas"] <- 3



# =================================================================================
# Summary Statistics
# =================================================================================

california <- subset(data,state=="California")
others <- subset(data,state!="California")
california.pre <- data[which(data$state=="California" & data$year <= 1988),]
california.post <- data[which(data$state=="California" & data$year > 1988),]
others.pre <- data[which(data$state!="California" & data$year <= 1988),]
others.post <- data[which(data$state!="California" & data$year > 1988),]

stat.desc(california.pre)
stat.desc(california.post)
stat.desc(others.pre)
stat.desc(others.post)


# =================================================================================
# Replication: Conventional SCM (Abadie et al., 2010)
# =================================================================================

dataprep.sc <- dataprep(
  foo = data,
  predictors = c("lnincome", "beer", "age15to24", "retprice"),
  predictors.op = "mean",
  time.predictors.prior = 1970:1988,
  special.predictors = list(
    list("cigsale", 1988 ,"mean"),
    list("cigsale", 1980, "mean"),
    list("cigsale", 1975, "mean")),
  dependent = "cigsale",
  unit.variable = "stateno",
  unit.names.variable = "state",
  time.variable = "year",
  treatment.identifier = 1,
  controls.identifier = c(2:39),
  time.optimize.ssr = 1970:1988,
  time.plot = 1970:2000)

synth.out.sc <- synth(data.prep.obj = dataprep.sc, method = "BFGS")

gaps <- dataprep.sc$Y1plot - (dataprep.sc$Y0plot %*% synth.out.sc$solution.w)

png(file = "Figure 1.png")
path.plot(synth.res = synth.out.sc, dataprep.res = dataprep.sc,
          Ylab = "CigSale", Xlab = "year",
          Ylim = c(0, 200), Legend = c("California",
                                      "synthetic California"), 
          Legend.position = "bottomright")
dev.off()

error.scm <- round(mean(gaps[1:18]^2), 4)


# ==================================================================================
# Tree-based Estimation
# =================================================================================

B <- 100

# Function of Tree-based Synthetic Control Method
tree.sc <- function(data, B, type = 1){
  
  output <- matrix(nrow = length(unique(data$year)), ncol = B+2)
  
  for (i in 1:B) {
    
    # Generate dataframes for trees
    set.seed(i)
    data.pre <- data[data$year <= 1988,]
    data.post <- dplyr::anti_join(data, data.pre)
    data.split <- initial_split(data.pre, prop = .67)
    data.train <- training(data.split)
    data.test <- testing(data.split)
    
    # Set the default regression to be fixed-effect regression, otherwise linear regression.
    if (type == 1) {
      formula.tree <- cigsale ~ lnincome + beer + age15to24 + retprice - 1
    } else {
      formula.tree <- cigsale ~ lnincome + beer + age15to24 + retprice
    }
    
    tree <- rpart(formula = formula.tree, data = data.train, method = "anova", control = list(maxdepth = 3))
    # Tune complexity parameter
    min_cp <- tree$cptable[which.min(tree$cptable[, "xerror"]), "CP"]
    tree.prune <- prune(tree, cp = min_cp)
    #rpart.plot(tree.prune)
    data.train$where <- tree.prune$where
    
    year.index <- max(data.train$year[data.train$stateno == 1])
    tree.no <- as.numeric(data.train$where[data.train$stateno == 1 & data.train$year == year.index])
    donor.pool <- data.train[data.train$where == tree.no,]
    data.post$where <- tree.no
    index <- c(unique(donor.pool$stateno))
    data.sc <- subset(data, stateno %in% index)
    
    dataprep.tree <- dataprep(
      foo = data.sc,
      predictors = c("lnincome", "beer", "age15to24", "retprice"),
      predictors.op = "mean",
      time.predictors.prior = 1970:1988,
      special.predictors = list(
        list("cigsale", 1988 ,"mean"),
        list("cigsale", 1980, "mean"),
        list("cigsale", 1975, "mean")),
      dependent = "cigsale",
      unit.variable = "stateno",
      unit.names.variable = "state",
      time.variable = "year",
      treatment.identifier = 1,
      controls.identifier = c(sort(unique(donor.pool$stateno)))[2:length(unique(donor.pool$stateno))],
      time.optimize.ssr = 1970:1988,
      time.plot = 1970:2000)
    
    
    synth.out.tree <- synth(data.prep.obj = dataprep.tree, method = "BFGS")
    output[,i] <- matrix(dataprep.tree$Y0plot %*% synth.out.tree$solution.w, ncol = 1)

  }
  
  output[,B+1] <- matrix(rowMins(output, value = TRUE), ncol = 1)
  output[,B+2] <- matrix(rowMaxs(output, value = TRUE), ncol = 1)
  output <- data.frame(output)
  output$year <- c(1970:2000)
  
  return(output)
  
}

tree.scm.lr <- tree.sc(data, B, 0)
tree.scm.fe <- tree.sc(data, B, 1)
# tree.scm.fe.alabama <- tree.sc(data, B, 1)

# Pre-treatment Prediciton Error
gaps.tree <- dataprep.sc$Y1plot - tree.scm.fe[,1:B]
gaps.tree.sq <- as.matrix(gaps.tree^2)
tt <- colMeans(gaps.tree.sq)
gaps.tree.mean <- mean(gaps.tree.sq[1:18,])
# Prediction Interval Width
interval.width <- matrix(tree.scm.fe[,B+2] - tree.scm.fe[,B+1], ncol = 1)

# ===========================================================================

# Figure 2
png(file = "Figure 2.png",width = 1000, height = 1000)
par(mfrow = c(2,2))

# Subplot 1
path.plot(synth.res = synth.out.sc, dataprep.res = dataprep.sc,
          Ylab = "CigSale", Xlab = "year",
          Ylim = c(0, 180),
          legend("topleft",
                 legend = c("California", "synthetic California", "Interval"),
                 lty = c(1,2,1),
                 col = c("black", "black", "grey"),
                 pt.cex = 1, 
                 cex = 0.6),
          Main = "Tree-based SCM (FE)")
abline(v = 1988, lty = 3, col = "black")

for (j in 1:B) {
  lines(tree.scm.fe$year, tree.scm.fe[,j], type = "l", lty = 3,  col = "grey")
}

# Subplot 2
path.plot(synth.res = synth.out.sc, dataprep.res = dataprep.sc,
          Ylab = "CigSale", Xlab = "year",
          Ylim = c(0, 180),
          legend("topleft",
                 legend = c("California", "synthetic California", "Interval"),
                 lty = c(1,2,1),
                 col = c("black", "black", "grey"),
                 pt.cex = 1, 
                 cex = 0.6),
          Main = "Tree-based SCM (OLS)")
abline(v = 1988, lty = 3, col = "black")

for (j in 1:B) {
  lines(tree.scm.lr$year, tree.scm.lr[,j], type = "l", lty = 3, col = "grey")
}

# Subplot 3 (The following codes are also used to plot Figure 3.)
path.plot(synth.res = synth.out.sc, dataprep.res = dataprep.sc,
          Ylab = "CigSale", Xlab = "year",
          Ylim = c(0, 180),
          legend("topleft",
                 legend = c("California", "synthetic California", "Interval"),
                #legend = c("Alabama", "synthetic Alabama", "Interval"),
                 lty = c(1,2,1),
                 col = c("black", "black", "grey"),
                 pt.cex = 1, 
                 cex = 0.6),
          Main = "Tree-based SCM (FE)")
abline(v = 1988, lty = 3, col = "black")

lines(tree.scm.fe$year, tree.scm.fe[,B+1], type = "l", lty = 3,  col = "red")
lines(tree.scm.fe$year, tree.scm.fe[,B+2], type = "l", lty = 3,  col = "red")
polygon(c(tree.scm.fe$year,rev(tree.scm.fe$year)),c(tree.scm.fe[,B+2],rev(tree.scm.fe[,B+1])),
        col = adjustcolor("grey",alpha.f = 0.5), border = NA)


# Subplot 4
path.plot(synth.res = synth.out.sc, dataprep.res = dataprep.sc,
          Ylab = "CigSale", Xlab = "year",
          Ylim = c(0, 180),
          legend("topleft",
                 legend = c("California", "synthetic California", "Interval"),
                 lty = c(1,2,1),
                 col = c("black", "black", "grey"),
                 pt.cex = 1, 
                 cex = 0.6),
          Main = "Tree-based SCM (OLS)")
abline(v = 1988, lty = 3, col = "black")

lines(tree.scm.lr$year, tree.scm.lr[,B+1], type = "l", lty = 3,  col = "blue")
lines(tree.scm.lr$year, tree.scm.lr[,B+2], type = "l", lty = 3,  col = "blue")
polygon(c(tree.scm.lr$year,rev(tree.scm.lr$year)),c(tree.scm.lr[,B+2],rev(tree.scm.lr[,B+1])),
        col = adjustcolor("grey",alpha.f = 0.5), border = NA)


dev.off()


# ================================================================================
# Other SC Competitors
# =================================================================================

# 1. Generalized SCM (Xu, 2017)
test < -regsynth(dataprep.sc$X0, dataprep.sc$X1, dataprep.sc$Y0plot[31,], dataprep.sc$Y1plot[31,], 1)

data.gscm <- data
data.gscm$D <- ifelse(data.gscm$stateno == 1 & data.gscm$year >= 1988, 1, 0)
g.scm <- gsynth(cigsale ~ D + lnincome, data = data.gscm, index = c("state","year"), 
              force = "two-way", CV = TRUE, r = c(0, 5), se = TRUE, 
              inference = "parametric", nboots = 1000, parallel = FALSE)

plot(tree.scm.lr$year, g.scm$Y, type = "l", lty = 1,  col = "black")
lines(tree.scm.fe$year, tree.scm.fe[,B+1], type = "l", lty = 3,  col = "red")
lines(tree.scm.fe$year, tree.scm.fe[,B+2], type = "l", lty = 3,  col = "red")
lines(tree.scm.fe$year, dataprep.sc$Y1plot, type = "l", lty = 3,  col = "red")


# 2. Augmented SCM (Ben-Michael et al., 2020)

aug.scm <- augsynth(cigsale ~ D,
                    stateno, year, t_int = 1988, data.gscm, progfunc = "None", scm = T)
aug.scm.pred <- matrix(c(1970:2000), ncol = 1)
aug.scm.pred <- cbind(aug.scm.pred, 
                      matrix(aug.scm[["data"]][["synth_data"]][["Y0plot"]] %*% aug.scm[["weights"]], ncol = 1))
scm.pred     <- matrix(dataprep.sc$Y0plot %*% synth.out.sc$solution.w, ncol = 1)

# Figure 4
png(file = "Figure 3.png", width = 1000, height = 1000)

plot(aug.scm.pred[,1], aug.scm.pred[,2], type = "l", lty = 2, xlab = "year",
     ylab = "Cigsale")
lines(tree.scm.fe$year, scm.pred, type = "l", lty = 2, col = "blue")
lines(tree.scm.fe$year, tree.scm.fe[,B+1], type = "l", lty = 3,  col = "red")
lines(tree.scm.fe$year, tree.scm.fe[,B+2], type = "l", lty = 3,  col = "red")
polygon(c(tree.scm.fe$year,rev(tree.scm.fe$year)),c(tree.scm.fe[,B+2],rev(tree.scm.fe[,B+1])),
        col = adjustcolor("grey",alpha.f = 0.5), border = NA)
abline(v = 1988, lty = 3, col = "black")

dev.off()



