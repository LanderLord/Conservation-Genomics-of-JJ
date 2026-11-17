#Notes:as we analyzed brown bear (Ursus arctos) populations from two regions (Himalaya and Tibet), for simplicity, this script will just use the 33 sampled populations as an example.
#This script used for Landscape connectivity analysis.
#The first part shows how to use ResistanceGA optimizing resistance surfaces, and the second part used for MLPE
############
#first part#
############
#install.packages("devtools")
#devtools::install_github("wpeterman/ResistanceGA") # v4.2-2
library(ResistanceGA)
library(JuliaCall)
library(raster)
library(sp)
library(doParallel)
#Before optimizing, you need to stack all environment layers (bio2, bio3, bio8, bio9, bio13, bio14, bio15, bio18, bio19, elev) into one ASCII folder.
#For more details, plz refer to https://github.com/wpeterman/ResistanceGA/
#Load sample coordinates and genetic distance
env <- read.csv("GEBV_WITH_LONLAT.csv", header = TRUE)
sp.dat <- SpatialPoints(env[, c("lon", "lat")])
geno <- t(read.table("417snps_genotype_matrix.tsv", sep = ",",
                     header = TRUE, row.names = 1, check.names = FALSE))
gen.dist <- as.matrix(dist(geno))
#Optimization with CIRCUITSCAPE (v4) in julia platform
JULIA_HOME <- "D:/Julia-1.7.0/bin"
julia_setup(JULIA_HOME)
#set threads number you want to use
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
jl.inputs <- jl.prep(n.Pops = 33,
                     response = as.vector(gen.dist[lower.tri(gen.dist)]),
                     CS_Point.File = sp.dat,
                     JULIA_HOME = JULIA_HOME,
                     run_test = TRUE)
GA.inputs <- GA.prep(ASCII.dir = "D:/connectivity/asc/",
                     Results.dir = "D:/connectivity/Result/",
                     min.cat = 1,
                     max.cat = 100,
                     max.cont = 100,
                     method = "AIC",
                     seed = 555,
                     parallel = cl,
                     run = 5,
                     maxiter = 10)
#Export info to cluster
clusterExport(cl = cl, varlist = c("jl.inputs", "GA.inputs", "sp.dat", "gen.dist"))
clusterCall(cl = cl, library, package = "ResistanceGA", character.only = TRUE)
clusterCall(cl = cl, library, package = "raster", character.only = TRUE)
jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)
stopCluster(cl)
#############
#second part#
#############
#install.packages("lme4")
library(lme4) #v1.1-27.1
library(Matrix)
library(geosphere)
#read sample info, genetic distance and geographic (Haversine) distance
info <- read.csv("GEBV_WITH_LONLAT.csv", stringsAsFactors = FALSE)[, c("SampleID", "lon", "lat", "population")]
rownames(info) <- info$SampleID
geno <- t(read.table("417snps_genotype_matrix.tsv", sep = ",",
                     header = TRUE, row.names = 1, check.names = FALSE))
geno <- geno[rownames(info), ]
gen.dist <- as.matrix(dist(geno))
fst <- gen.dist[lower.tri(gen.dist)]
geo <- distm(info[, c("lon", "lat")], fun = distHaversine)
dis <- geo[lower.tri(geo)]
#Read resistance matrices optimized above
vars <- c("bio2", "bio3", "bio8", "bio9", "bio13", "bio14", "bio15", "bio18", "bio19", "elev")
res_list <- lapply(vars, function(x){
  m <- as.matrix(read.csv(paste0("./Results/", x, "_jlResistMat.csv"), header = FALSE))
  m[lower.tri(m)]
})
names(res_list) <- vars
#Keep Himalaya-Tibet cross-region pairs only
pair <- which(lower.tri(matrix(TRUE, nrow(info), nrow(info))), arr.ind = TRUE)
pop1 <- info$population[pair[, 1]]
pop2 <- info$population[pair[, 2]]
cross <- (pop1 == "Himalaya" & pop2 == "Tibet") | (pop1 == "Tibet" & pop2 == "Himalaya")
df <- data.frame(Genetic_distance = fst, res_list, dis = dis, pop1 = pop1, pop2 = pop2)[cross, ]
#standardization
df[, c(vars, "dis")] <- lapply(df[, c(vars, "dis")], scale)
#Build the sparse random-effect matrix for MLPE (population as random effect)
lev <- sort(unique(c(df$pop1, df$pop2)))
ZZ <- Reduce("+", lapply(list(df$pop1, df$pop2), function(x){
  Matrix::fac2sparse(factor(x, levels = lev), "d", drop = FALSE)
}))
df$MLPE_ID <- factor(rep(lev, length.out = nrow(df)), levels = lev)
#Define the MLPE function (REML = FALSE for AIC comparison)
MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  lme4::mkMerMod(environment(dfun), opt, mod2$reTrms, fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  #refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  lme4::mkMerMod(environment(dfun), opt, mod2$reTrms, fr = mod2$fr)
}
#Models: geographic, single bio, temperature, humidity, annual variability, elevation, climate, climate + elevation
forms <- list(
  Geographic_distance = Genetic_distance ~ dis + (1|MLPE_ID),
  bio2 = Genetic_distance ~ bio2 + (1|MLPE_ID),
  bio3 = Genetic_distance ~ bio3 + (1|MLPE_ID),
  bio8 = Genetic_distance ~ bio8 + (1|MLPE_ID),
  bio9 = Genetic_distance ~ bio9 + (1|MLPE_ID),
  bio13 = Genetic_distance ~ bio13 + (1|MLPE_ID),
  bio14 = Genetic_distance ~ bio14 + (1|MLPE_ID),
  bio15 = Genetic_distance ~ bio15 + (1|MLPE_ID),
  bio18 = Genetic_distance ~ bio18 + (1|MLPE_ID),
  bio19 = Genetic_distance ~ bio19 + (1|MLPE_ID),
  Elevation = Genetic_distance ~ elev + (1|MLPE_ID),
  Temperature = Genetic_distance ~ bio2 + bio3 + bio8 + bio9 + (1|MLPE_ID),
  Precipitation = Genetic_distance ~ bio13 + bio14 + bio15 + bio18 + bio19 + (1|MLPE_ID),
  Annual_variability = Genetic_distance ~ bio2 + bio3 + bio15 + (1|MLPE_ID),
  Climate = Genetic_distance ~ bio2 + bio3 + bio8 + bio9 + bio13 + bio14 + bio15 + bio18 + bio19 + (1|MLPE_ID),
  Climate_Elevation = Genetic_distance ~ bio2 + bio3 + bio8 + bio9 + bio13 + bio14 + bio15 + bio18 + bio19 + elev + (1|MLPE_ID)
)
mods <- lapply(forms, MLPEnoREML, data = df)
#Compute AICc and the marginal/conditional R2
res <- data.frame(Model = names(mods),
                  K = sapply(mods, function(m) attr(logLik(m), "df")),
                  AICc = sapply(mods, function(m){
                    k <- attr(logLik(m), "df"); n <- nobs(m)
                    AIC(m) + 2 * k * (k + 1) / (n - k - 1)
                  }),
                  row.names = NULL)
res$Delta_AICc <- res$AICc - min(res$AICc)
res$Akaike_weight <- exp(-0.5 * res$Delta_AICc) / sum(exp(-0.5 * res$Delta_AICc))
res <- res[order(res$AICc), ]
#Save results
write.table(res, "MLPE_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
res[1, ]
