# =============================================================================
# SPATIAL ECONOMETRICS PROJECT – Rise of Chega
# =============================================================================

# Author: Anil Mujagic
# Student-ID: 65067


#==============================================================================#
# Packages                                                                     #
#==============================================================================#

# install.packages(c("spdep","sphet","spatialreg","psych","MASS","sf",
#                    "mapview","lmtest","spgwr"))

library(spdep)
library(sphet)
library(spatialreg)
library(psych)
library(MASS)
library(sf)
library(mapview)
library(lmtest)
library(spgwr)

rm(list = ls())

setwd("C:\\Users\\anilm\\Documents\\Chega Spatial\\Data")


# =============================================================================
# HELPERS – fault-tolerant wrappers
# =============================================================================

try_run <- function(label, expr) {
  tryCatch(expr,
           error = function(e) {
             message("SKIPPED [", label, "]: ", conditionMessage(e))
             NULL
           })
}

safe_summary <- function(x, label = "") {
  if (is.null(x)) {
    message("SKIPPED summary [", label, "]: model not available")
    return(invisible(NULL))
  }
  tryCatch(summary(x),
           error = function(e) {
             message("SKIPPED summary [", label, "]: ", conditionMessage(e))
             invisible(NULL)
           })
}

safe_impacts <- function(model, label, listw, R = 1000) {
  if (is.null(model)) {
    message("SKIPPED impacts [", label, "]: model not available")
    return(NULL)
  }
  tryCatch(
    summary(impacts(model, listw = listw, R = R), zstats = TRUE),
    error = function(e) {
      message("SKIPPED impacts [", label, "]: ", conditionMessage(e))
      NULL
    }
  )
}

pack_peffs <- function(peffs) {
  if (is.null(peffs)) return(NULL)
  out <- list(
    effects  = `colnames<-`(cbind(peffs$res$direct,
                                  peffs$res$indirect,
                                  peffs$res$total),
                            c("Direct", "Indirect", "Total")),
    std.err  = peffs$semat,
    p.values = round(peffs$pzmat, 4)
  )
  rownames(out$effects) <- rownames(peffs$semat)
  out
}

safe_lr <- function(m1, m2, label) {
  if (is.null(m1) || is.null(m2)) {
    message("SKIPPED LR test [", label, "]: one or both models not available")
    return(invisible(NULL))
  }
  tryCatch({ cat("\n--- LR test:", label, "---\n"); print(lrtest(m1, m2)) },
           error = function(e) message("SKIPPED LR [", label, "]: ", conditionMessage(e)))
}

safe_moran_resid <- function(model, label, listw) {
  if (is.null(model)) {
    message("SKIPPED Moran resid [", label, "]: model not available")
    return(NULL)
  }
  tryCatch({
    resid <- as.numeric(residuals(model))
    mt    <- moran.test(resid, listw = listw, zero.policy = TRUE)
    cat(sprintf("  %-38s  I = %6.4f  z = %6.3f  p = %.4f\n",
                label, mt$estimate[1], mt$statistic, mt$p.value))
    mt
  }, error = function(e) {
    message("SKIPPED Moran resid [", label, "]: ", conditionMessage(e))
    NULL
  })
}

sphet_pv <- function(m) {
  tryCatch(
    as.numeric(2 * (1 - pnorm(abs(m$coef / sqrt(diag(m$var)))))),
    error = function(e) rep(NA_real_, length(m$coef))
  )
}

# coef_vec(): force sphet $coef to a plain named numeric vector.
# sphet occasionally returns a list or matrix instead of a vector,
# which causes "invalid type list in rbind.Matrix" errors.
coef_vec <- function(m) as.numeric(unlist(m$coef))

na_block <- function(nrow) matrix(NA_real_, nrow, 2)

safe_block <- function(nrow, label, expr) {
  tryCatch(expr,
           error = function(e) {
             message("SKIPPED comp block [", label, "]: ", conditionMessage(e))
             na_block(nrow)
           })
}

safe_bp <- function(model, label, sarlm = FALSE) {
  if (is.null(model)) {
    message("SKIPPED BP [", label, "]: model not available"); return(NULL)
  }
  tryCatch(
    if (sarlm) bptest.Sarlm(model, studentize = TRUE)
    else       bptest(model,       studentize = TRUE),
    error = function(e) { message("SKIPPED BP [", label, "]: ", conditionMessage(e)); NULL }
  )
}

safe_anselin <- function(bp, lm_comp, label) {
  if (is.null(bp) || is.null(lm_comp)) {
    message("SKIPPED Anselin [", label, "]: BP or LM test not available")
    return(invisible(NULL))
  }
  tryCatch({
    stat <- bp$statistic + lm_comp$statistic
    df   <- bp$parameter + lm_comp$parameter
    out  <- c(statistic = stat, df = df, `p-value` = 1 - pchisq(stat, df))
    cat("\n--- Anselin joint test:", label, "---\n")
    print(round(out, 4))
  }, error = function(e) message("SKIPPED Anselin [", label, "]: ", conditionMessage(e)))
}


# =============================================================================
# DATA LOADING
# =============================================================================
# politics.shp: CHC (dependent) | Utopia, Ghost.Town, Tik.T_.Gen, Chega_aign

data <- st_read("politics.shp", quiet = TRUE)


# =============================================================================
# POLYGONS & CENTROIDS
# =============================================================================

poly      <- st_geometry(data)
centroids <- st_centroid(poly)
coords    <- st_coordinates(centroids)

plot(poly, col = "lightblue", border = "black", main = "Polygons with Centroids")
plot(centroids, col = "red", pch = 20, cex = 1.2, add = TRUE)


# =============================================================================
# SPATIAL WEIGHTS MATRICES
# =============================================================================

## Distance-band (10 km)
nn <- dnearneigh(coords, 0, 10, longlat = TRUE)
plot(poly, col = "lightblue", border = "black",
     main = "Distance-band neighbours (10 km)")
plot(centroids, col = "red", pch = 20, cex = 1.2, add = TRUE)
plot(nn, coords, col = "blue", add = TRUE)

## k-Nearest Neighbour (k = 5)
knn     <- knearneigh(coords, 5, longlat = TRUE)
knn_nbl <- knn2nb(knn)
plot(poly, col = "lightblue", border = "black", main = "5-NN neighbours")
plot(centroids, col = "red", pch = 20, cex = 1.2, add = TRUE)
plot(knn_nbl, coords, col = "blue", add = TRUE)

## Queen contiguity [used for all estimation]
queen <- poly2nb(poly)
plot(poly, col = "lightblue", border = "black",
     main = "Queen contiguity neighbours")
plot(centroids, col = "red", pch = 20, cex = 1.2, add = TRUE)
plot(queen, coords, col = "blue", add = TRUE)

## Rook contiguity
rook <- poly2nb(poly, queen = FALSE)
plot(poly, col = "lightblue", border = "black",
     main = "Rook contiguity neighbours")
plot(centroids, col = "red", pch = 20, cex = 1.2, add = TRUE)
plot(rook, coords, col = "blue", add = TRUE)

wmat <- nb2listw(queen, style = "W")


# =============================================================================
# GLOBAL & LOCAL MORAN'S I (dependent variable)
# =============================================================================

## Global Moran's I – randomisation
gmoran <- moran.test(data$CHC, wmat, alternative = "two.sided")
print(gmoran)

## Global Moran's I – normality assumption
gmoran_norm <- moran.test(data$CHC, wmat, alternative = "two.sided",
                          randomisation = FALSE)
print(gmoran_norm)

## Local Moran's I
lmoran    <- localmoran(data$CHC, wmat, alternative = "two.sided")
lmoran_df <- as.data.frame(lmoran)
lmoran_df <- lmoran_df[order(-lmoran_df$Ii), ]
hist(lmoran_df$Ii, main = "Local Moran's I – CHC")

## Check mean of local = global
gmoran$estimate[1] == mean(lmoran_df$Ii)

## Global Moran plot
moran.plot(data$CHC, wmat)

## Local Moran plot
demean_CHC <- data$CHC - mean(data$CHC)
sig        <- 0.1

norm_fact_lmoran <- (nrow(lmoran) / sum(nb2mat(queen, style = "W")))
mean_lmoran      <- norm_fact_lmoran * mean(lmoran[, 1])
demean_lmoran    <- lmoran[, 1] - mean_lmoran

quadrant <- rep(0, nrow(data))
quadrant[demean_CHC > 0 & demean_lmoran > 0] <- 4   # High-High
quadrant[demean_CHC < 0 & demean_lmoran < 0] <- 1   # Low-Low
quadrant[demean_CHC < 0 & demean_lmoran > 0] <- 3   # High-Low
quadrant[demean_CHC > 0 & demean_lmoran < 0] <- 2   # Low-High
quadrant[lmoran[, 5] > sig]                  <- 0   # Non-significant

brks    <- seq(0, 4, 1)
colors1 <- c("white", "blue", rgb(0,0,1,alpha=0.4), rgb(1,0,0,alpha=0.4), "red")
par(oma = c(0, 0, 3, 5))
plot(poly, col = colors1[findInterval(quadrant, brks, all.inside = FALSE)],
     main = "Local Moran's I – CHC")
legend("right",
       legend = c("insignificant","low-low","low-high","high-low","high-high"),
       fill = colors1, bty = "n", cex = 0.8)


# =============================================================================
# VARIABLE SETUP & SPATIAL LAGS
# =============================================================================

names_x  <- c("Utopia", "Ghost.Town", "Tik.T_.Gen", "Chega_aign")

## Generate spatially lagged X variables (same approach as professor)
for (v in names_x) {
  data[[paste0("W_", v)]] <- lag.listw(wmat, as.numeric(data[[v]]),
                                       zero.policy = TRUE)
}

names_wx <- paste0("W_", names_x)

## Formulas
### spec    : baseline — only own-unit regressors X
### spec_wx : Durbin  — X plus spatially-lagged WX
###
### IMPORTANT: Durbin models are estimated by passing spec_wx directly.
### Do NOT use Durbin = TRUE inside spreg/errorsarlm/lagsarlm/sacsarlm:
### it produces a non-symmetric variance matrix that breaks impacts()
### with "sigma must be a symmetric matrix".

spec    <- as.formula(paste("CHC ~", paste(names_x,  collapse = " + ")))
spec_wx <- as.formula(paste("CHC ~", paste(names_x,  collapse = " + "),
                             "+",    paste(names_wx, collapse = " + ")))


# =============================================================================
# ESTIMATION
# =============================================================================

# ---- OLS --------------------------------------------------------------------
ols    <- try_run("OLS",    spreg(spec,    listw = wmat, model = "ols", data = data))
ols_wx <- try_run("OLS WX", spreg(spec_wx, listw = wmat, model = "ols", data = data))
safe_summary(ols,    "OLS")
safe_summary(ols_wx, "OLS WX")

## lm objects needed for AIC/BIC/LR/BP/RESET tests
ols_lm    <- try_run("OLS lm",    lm(spec,    data = data))
ols_wx_lm <- try_run("OLS WX lm", lm(spec_wx, data = data))


# ---- SEM – Maximum Likelihood -----------------------------------------------
sem_ml    <- try_run("SEM ML",    errorsarlm(spec,    listw = wmat, zero.policy = TRUE, data = data))
sem_wx_ml <- try_run("SEM WX ML", errorsarlm(spec_wx, listw = wmat, zero.policy = TRUE, data = data))
safe_summary(sem_ml,    "SEM ML")
safe_summary(sem_wx_ml, "SEM WX ML")

sem_wx_ml_peffs    <- safe_impacts(sem_wx_ml, "SEM WX ML", wmat)
sem_wx_ml_peffs_pv <- pack_peffs(sem_wx_ml_peffs)


# ---- SEM – KP FGLS ----------------------------------------------------------
#### Note: "rho" in spreg output = lambda (spatial error parameter)
sem_gls    <- try_run("SEM GLS",    spreg(spec,    listw = wmat, model = "error", data = data))
sem_wx_gls <- try_run("SEM WX GLS", spreg(spec_wx, listw = wmat, model = "error", data = data))
safe_summary(sem_gls,    "SEM GLS")
safe_summary(sem_wx_gls, "SEM WX GLS")

sem_wx_gls_peffs    <- safe_impacts(sem_wx_gls, "SEM WX GLS", wmat)
sem_wx_gls_peffs_pv <- pack_peffs(sem_wx_gls_peffs)


# ---- SEM – KP Robust [PREFERRED] -------------------------------------------
#### het = TRUE: Kelejian-Prucha (2010) heteroskedasticity-robust variance
sem_kp    <- try_run("SEM KP",    spreg(spec,    listw = wmat, model = "error", het = TRUE, data = data))
sem_wx_kp <- try_run("SEM WX KP", spreg(spec_wx, listw = wmat, model = "error", het = TRUE, data = data))
safe_summary(sem_kp,    "SEM KP")
safe_summary(sem_wx_kp, "SEM WX KP")

sem_wx_kp_peffs    <- safe_impacts(sem_wx_kp, "SEM WX KP", wmat)
sem_wx_kp_peffs_pv <- pack_peffs(sem_wx_kp_peffs)


# ---- SLM – Maximum Likelihood -----------------------------------------------
slm_ml    <- try_run("SLM ML",    lagsarlm(spec,    listw = wmat, data = data))
slm_wx_ml <- try_run("SLM WX ML", lagsarlm(spec_wx, listw = wmat, data = data))
safe_summary(slm_ml,    "SLM ML")
safe_summary(slm_wx_ml, "SLM WX ML")

slm_ml_peffs       <- safe_impacts(slm_ml,    "SLM ML",    wmat)
slm_ml_peffs_pv    <- pack_peffs(slm_ml_peffs)
slm_wx_ml_peffs    <- safe_impacts(slm_wx_ml, "SLM WX ML", wmat)
slm_wx_ml_peffs_pv <- pack_peffs(slm_wx_ml_peffs)


# ---- SLM – Two-Stage Least Squares ------------------------------------------
#### Note: "lambda" in spreg output = rho (spatial lag parameter)
slm_2sls    <- try_run("SLM 2SLS",    spreg(spec,    listw = wmat, model = "lag", data = data))
slm_wx_2sls <- try_run("SLM WX 2SLS", spreg(spec_wx, listw = wmat, model = "lag", data = data))
safe_summary(slm_2sls,    "SLM 2SLS")
safe_summary(slm_wx_2sls, "SLM WX 2SLS")

slm_2sls_peffs       <- safe_impacts(slm_2sls,    "SLM 2SLS",    wmat)
slm_2sls_peffs_pv    <- pack_peffs(slm_2sls_peffs)
slm_wx_2sls_peffs    <- safe_impacts(slm_wx_2sls, "SLM WX 2SLS", wmat)
slm_wx_2sls_peffs_pv <- pack_peffs(slm_wx_2sls_peffs)


# ---- SLM – White robust -----------------------------------------------------
slm_white    <- try_run("SLM White",    spreg(spec,    listw = wmat, model = "lag", het = TRUE, data = data))
slm_wx_white <- try_run("SLM WX White", spreg(spec_wx, listw = wmat, model = "lag", het = TRUE, data = data))
safe_summary(slm_white,    "SLM White")
safe_summary(slm_wx_white, "SLM WX White")

slm_white_peffs       <- safe_impacts(slm_white,    "SLM White",    wmat)
slm_white_peffs_pv    <- pack_peffs(slm_white_peffs)
slm_wx_white_peffs    <- safe_impacts(slm_wx_white, "SLM WX White", wmat)
slm_wx_white_peffs_pv <- pack_peffs(slm_wx_white_peffs)


# ---- SARAR – Maximum Likelihood ---------------------------------------------
sarar_ml    <- try_run("SARAR ML",    sacsarlm(spec,    listw = wmat, data = data))
sarar_wx_ml <- try_run("SARAR WX ML", sacsarlm(spec_wx, listw = wmat, data = data))
safe_summary(sarar_ml,    "SARAR ML")
safe_summary(sarar_wx_ml, "SARAR WX ML")

sarar_ml_peffs       <- safe_impacts(sarar_ml,    "SARAR ML",    wmat)
sarar_ml_peffs_pv    <- pack_peffs(sarar_ml_peffs)
sarar_wx_ml_peffs    <- safe_impacts(sarar_wx_ml, "SARAR WX ML", wmat)
sarar_wx_ml_peffs_pv <- pack_peffs(sarar_wx_ml_peffs)


# ---- OLS White --------------------------------------------------------------
ols_white    <- try_run("OLS White",    spreg(spec,    listw = wmat, model = "ols", het = TRUE, data = data))
ols_wx_white <- try_run("OLS WX White", spreg(spec_wx, listw = wmat, model = "ols", het = TRUE, data = data))
safe_summary(ols_white,    "OLS White")
safe_summary(ols_wx_white, "OLS WX White")


# =============================================================================
# COEFFICIENT COMPARISON TABLES
# =============================================================================

## --- No spatial lags of X ---------------------------------------------------

ols_comp <- safe_block(6, "OLS", {
  if (!is.null(ols)) rbind(na_block(2), cbind(coef_vec(ols), sphet_pv(ols)))
  else na_block(6)
})
colnames(ols_comp) <- c("OLS", "p-value")

sem_ml_comp <- safe_block(6, "SEM ML", {
  if (!is.null(sem_ml))
    rbind(na_block(1),
          cbind(c(sem_ml$lambda, sem_ml$coefficients),
                c(2*(1-pnorm(abs(sem_ml$lambda / sem_ml$lambda.se))),
                  2*(1-pnorm(abs(sem_ml$coefficients / sem_ml$rest.se))))))
  else na_block(6)
})
colnames(sem_ml_comp) <- c("SEM ML", "p-value")

sem_gls_comp <- safe_block(6, "SEM GLS", {
  if (!is.null(sem_gls)) {
    cv <- coef_vec(sem_gls)
    pv <- sphet_pv(sem_gls)
    rbind(na_block(1),
          cbind(c(tail(cv, 1), head(cv, -1)),
                c(tail(pv, 1), head(pv, -1))))
  } else na_block(6)
})
colnames(sem_gls_comp) <- c("SEM GLS", "p-value")

slm_ml_comp <- safe_block(6, "SLM ML", {
  if (!is.null(slm_ml))
    rbind(c(slm_ml$rho, 2*(1-pnorm(abs(slm_ml$rho / slm_ml$rho.se)))),
          na_block(1),
          cbind(slm_ml$coefficients,
                2*(1-pnorm(abs(slm_ml$coefficients / slm_ml$rest.se)))))
  else na_block(6)
})
colnames(slm_ml_comp) <- c("Spatial Lag ML", "p-value")

slm_2sls_comp <- safe_block(6, "SLM 2SLS", {
  if (!is.null(slm_2sls)) {
    cv <- coef_vec(slm_2sls)
    pv <- sphet_pv(slm_2sls)
    rbind(c(tail(cv, 1), tail(pv, 1)),
          na_block(1),
          cbind(head(cv, -1), head(pv, -1)))
  } else na_block(6)
})
colnames(slm_2sls_comp) <- c("Spatial Lag 2SLS", "p-value")

sarar_ml_comp <- safe_block(6, "SARAR ML", {
  if (!is.null(sarar_ml))
    rbind(c(sarar_ml$rho,    2*(1-pnorm(abs(sarar_ml$rho    / sarar_ml$rho.se)))),
          c(sarar_ml$lambda, 2*(1-pnorm(abs(sarar_ml$lambda / sarar_ml$lambda.se)))),
          cbind(sarar_ml$coefficients,
                2*(1-pnorm(abs(sarar_ml$coefficients / sarar_ml$rest.se)))))
  else na_block(6)
})
colnames(sarar_ml_comp) <- c("SARAR ML", "p-value")

comp <- tryCatch({
  out <- cbind(ols_comp, sem_ml_comp, sem_gls_comp,
               slm_ml_comp, slm_2sls_comp, sarar_ml_comp)
  out <- round(out, 3)
  rownames(out)[1:2] <- c("Rho", "Lambda")
  out
}, error = function(e) { message("SKIPPED comp table (no WX): ", conditionMessage(e)); NULL })

cat("\n--- Coefficient comparison (no WX) ---\n")
if (!is.null(comp)) print(comp)


## --- With spatial lags of X (Durbin) ----------------------------------------

ols_wx_comp <- safe_block(10, "OLS WX", {
  if (!is.null(ols_wx)) rbind(na_block(2), cbind(coef_vec(ols_wx), sphet_pv(ols_wx)))
  else na_block(10)
})
colnames(ols_wx_comp) <- c("OLS", "p-value")

sem_wx_ml_comp <- safe_block(10, "SEM WX ML", {
  if (!is.null(sem_wx_ml))
    rbind(na_block(1),
          cbind(c(sem_wx_ml$lambda, sem_wx_ml$coefficients),
                c(2*(1-pnorm(abs(sem_wx_ml$lambda / sem_wx_ml$lambda.se))),
                  2*(1-pnorm(abs(sem_wx_ml$coefficients / sem_wx_ml$rest.se))))))
  else na_block(10)
})
colnames(sem_wx_ml_comp) <- c("SEM ML", "p-value")

sem_wx_gls_comp <- safe_block(10, "SEM WX GLS", {
  if (!is.null(sem_wx_gls)) {
    cv <- coef_vec(sem_wx_gls)
    pv <- sphet_pv(sem_wx_gls)
    rbind(na_block(1),
          cbind(c(tail(cv, 1), head(cv, -1)),
                c(tail(pv, 1), head(pv, -1))))
  } else na_block(10)
})
colnames(sem_wx_gls_comp) <- c("SEM GLS", "p-value")

slm_wx_ml_comp <- safe_block(10, "SLM WX ML", {
  if (!is.null(slm_wx_ml))
    rbind(c(slm_wx_ml$rho, 2*(1-pnorm(abs(slm_wx_ml$rho / slm_wx_ml$rho.se)))),
          na_block(1),
          cbind(slm_wx_ml$coefficients,
                2*(1-pnorm(abs(slm_wx_ml$coefficients / slm_wx_ml$rest.se)))))
  else na_block(10)
})
colnames(slm_wx_ml_comp) <- c("SAR(1) ML", "p-value")

slm_wx_2sls_comp <- safe_block(10, "SLM WX 2SLS", {
  if (!is.null(slm_wx_2sls)) {
    cv <- coef_vec(slm_wx_2sls)
    pv <- sphet_pv(slm_wx_2sls)
    rbind(c(tail(cv, 1), tail(pv, 1)),
          na_block(1),
          cbind(head(cv, -1), head(pv, -1)))
  } else na_block(10)
})
colnames(slm_wx_2sls_comp) <- c("SAR(1) 2SLS", "p-value")

sarar_wx_ml_comp <- safe_block(10, "SARAR WX ML", {
  if (!is.null(sarar_wx_ml))
    rbind(c(sarar_wx_ml$rho,    2*(1-pnorm(abs(sarar_wx_ml$rho    / sarar_wx_ml$rho.se)))),
          c(sarar_wx_ml$lambda, 2*(1-pnorm(abs(sarar_wx_ml$lambda / sarar_wx_ml$lambda.se)))),
          cbind(sarar_wx_ml$coefficients,
                2*(1-pnorm(abs(sarar_wx_ml$coefficients / sarar_wx_ml$rest.se)))))
  else na_block(10)
})
colnames(sarar_wx_ml_comp) <- c("SARAR ML", "p-value")

comp_wx <- tryCatch({
  out <- cbind(ols_wx_comp, sem_wx_ml_comp, sem_wx_gls_comp,
               slm_wx_ml_comp, slm_wx_2sls_comp, sarar_wx_ml_comp)
  out <- round(out, 3)
  rownames(out)[1:2] <- c("Rho", "Lambda")
  out
}, error = function(e) { message("SKIPPED comp table (WX): ", conditionMessage(e)); NULL })

cat("\n--- Coefficient comparison (with WX) ---\n")
if (!is.null(comp_wx)) print(comp_wx)


# =============================================================================
# LM SPECIFICATION TESTS
# =============================================================================
# Run on BOTH specifications so we can compare how the spatial signal
# changes when WX is included — key for the revised model selection.

lm_tests <- if (!is.null(ols_lm)) {
  try_run("LM tests", lm.RStests(ols_lm, listw = wmat, test = "all",
                                 zero.policy = TRUE))
} else NULL

lm_tests_wx <- if (!is.null(ols_wx_lm)) {
  try_run("LM tests WX", lm.RStests(ols_wx_lm, listw = wmat, test = "all",
                                    zero.policy = TRUE))
} else NULL

if (!is.null(lm_tests))    { cat("\n--- LM tests (no WX) ---\n");   summary(lm_tests) }
if (!is.null(lm_tests_wx)) { cat("\n--- LM tests (with WX) ---\n"); summary(lm_tests_wx) }


# =============================================================================
# RAMSEY RESET TEST
# =============================================================================
# Tests linear functional form. Rejection would suggest non-linearities
# that could be masquerading as spatial dependence.

reset_ols    <- if (!is.null(ols_lm))    try_run("RESET",    resettest(ols_lm,    power = 2:3, type = "fitted")) else NULL
reset_ols_wx <- if (!is.null(ols_wx_lm)) try_run("RESET WX", resettest(ols_wx_lm, power = 2:3, type = "fitted")) else NULL

cat("\n--- Ramsey RESET test (no WX) ---\n")
if (!is.null(reset_ols))    print(reset_ols)    else message("SKIPPED")
cat("\n--- Ramsey RESET test (with WX) ---\n")
if (!is.null(reset_ols_wx)) print(reset_ols_wx) else message("SKIPPED")


# =============================================================================
# HETEROSKEDASTICITY TESTS (Spatial Breusch-Pagan)
# =============================================================================
# bptest()       on lm objects  (OLS)
# bptest.Sarlm() on Sarlm objects (all spatial ML models) = spatial BP test

bptest_ols      <- safe_bp(ols_lm,    "OLS")
bptest_ols_wx   <- safe_bp(ols_wx_lm, "OLS WX")
bptest_sem      <- safe_bp(sem_ml,    "SEM ML",     sarlm = TRUE)
bptest_sem_wx   <- safe_bp(sem_wx_ml, "SEM WX ML",  sarlm = TRUE)
bptest_slm      <- safe_bp(slm_ml,    "SLM ML",     sarlm = TRUE)
bptest_slm_wx   <- safe_bp(slm_wx_ml, "SLM WX ML",  sarlm = TRUE)
bptest_sarar    <- safe_bp(sarar_ml,  "SARAR ML",   sarlm = TRUE)
bptest_sarar_wx <- safe_bp(sarar_wx_ml,"SARAR WX ML",sarlm = TRUE)

cat("\n--- BP tests ---\n")
for (obj in list(bptest_ols, bptest_ols_wx,
                 bptest_sem, bptest_sem_wx,
                 bptest_slm, bptest_slm_wx,
                 bptest_sarar, bptest_sarar_wx)) {
  if (!is.null(obj)) print(obj)
}


# =============================================================================
# ANSELIN (1996) JOINT TEST
# =============================================================================

safe_anselin(bptest_sem,    lm_tests$RSerr,    "SEM (no WX)")
safe_anselin(bptest_slm,    lm_tests$RSlag,    "SLM (no WX)")
safe_anselin(bptest_sem_wx, lm_tests_wx$RSerr, "SEM (with WX)")
safe_anselin(bptest_slm_wx, lm_tests_wx$RSlag, "SLM (with WX)")


# =============================================================================
# GEOGRAPHICALLY WEIGHTED REGRESSION (GWR)
# =============================================================================

tryCatch({
  bandwidth <- gwr.sel(formula = spec, coords = coords, longlat = FALSE,
                       gweight = gwr.Gauss, adapt = TRUE, tol = 1e-06, data = data)
  gwr_model <- gwr(spec, coords = coords, adapt = bandwidth,
                   hatmatrix = TRUE, longlat = FALSE, data = data)
  gwr_df <- as(gwr_model$SDF, "data.frame")
  cat("\nGWR SDF column names:\n"); print(names(gwr_df))

  for (gwr_var in names_x) {
    estim_coef <- gwr_df[[gwr_var]]
    se_col     <- paste0(gwr_var, "_se")
    if (!se_col %in% names(gwr_df)) {
      warning(paste("GWR SE column not found:", se_col, "- skipping")); next
    }
    se_coef  <- gwr_df[[se_col]]
    pv_coef  <- 2 * (1 - pnorm(abs(estim_coef / se_coef)))
    sig_coef <- ifelse(pv_coef > 0.05, NA, estim_coef)

    describe(sig_coef)
    boxplot(sig_coef, main = paste("GWR coefficients:", gwr_var))
    hist(sig_coef,    main = paste("GWR coefficients:", gwr_var))

    brks_gwr <- c(-3,-2,-1,0,1,2,3) * sd(estim_coef, na.rm = TRUE)
    colors   <- rainbow(6)
    par(mar = c(0.1, 0.1, 4, 4))
    plot(poly, col = colors[findInterval(sig_coef, brks_gwr, all.inside = TRUE)],
         main = paste("GWR effects:", gwr_var), cex.main = 0.8)
    legend("right",
           legend = c("high<0","mod<0","low<0","low>0","mod>0","high>0"),
           fill = colors, bty = "n", cex = 0.8)
  }
}, error = function(e) message("SKIPPED GWR: ", conditionMessage(e)))


# =============================================================================
# FULL MODEL COMPARISON FRAMEWORK
# =============================================================================
# Addresses professor's feedback: test all candidate models on equal footing.
# The key insight is that LM tests on a baseline OLS without WX or Wy push
# all spatial signal into the error term by construction. The revised approach
# uses post-estimation Moran's I on residuals as the primary selection tool —
# the model whose residuals are spatially random has adequately captured the
# spatial structure of the data.


# --- 1. POST-ESTIMATION MORAN'S I ON RESIDUALS --------------------------------
# This is the primary model selection criterion.
# A model with significant residual Moran's I (p < 0.05) has NOT adequately
# captured the spatial structure, regardless of in-sample fit.

cat("\n")
cat("=============================================================================\n")
cat("FULL MODEL COMPARISON FRAMEWORK\n")
cat("=============================================================================\n")

cat("\n=== 1. POST-ESTIMATION MORAN'S I ON RESIDUALS ===\n")
cat("p < 0.05 = residual spatial autocorrelation remains (model inadequate)\n")
cat("p >= 0.05 = residuals spatially random (spatial structure captured)\n\n")
cat(sprintf("  %-38s  %s\n", "Model", "I        z        p-value"))
cat(paste(rep("-", 72), collapse = ""), "\n")

cat("\n  -- Baseline (no WX) --\n")
mr_ols         <- safe_moran_resid(ols_lm,     "OLS",                    wmat)
mr_sem_ml      <- safe_moran_resid(sem_ml,     "SEM ML",                 wmat)
mr_sem_gls     <- safe_moran_resid(sem_gls,    "SEM KP (FGLS)",          wmat)
mr_sem_kp      <- safe_moran_resid(sem_kp,     "SEM KP robust [pref]",   wmat)
mr_slm_ml      <- safe_moran_resid(slm_ml,     "SLM ML",                 wmat)
mr_slm_2sls    <- safe_moran_resid(slm_2sls,   "SLM 2SLS",               wmat)
mr_slm_white   <- safe_moran_resid(slm_white,  "SLM White",              wmat)
mr_sarar_ml    <- safe_moran_resid(sarar_ml,   "SARAR ML",               wmat)

cat("\n  -- Durbin / with WX --\n")
mr_ols_wx      <- safe_moran_resid(ols_wx_lm,    "OLS Durbin (SLX)",       wmat)
mr_sem_wx_ml   <- safe_moran_resid(sem_wx_ml,    "SEM ML Durbin",          wmat)
mr_sem_wx_gls  <- safe_moran_resid(sem_wx_gls,   "SEM KP Durbin (FGLS)",   wmat)
mr_sem_wx_kp   <- safe_moran_resid(sem_wx_kp,    "SEM KP robust Durbin",   wmat)
mr_slm_wx_ml   <- safe_moran_resid(slm_wx_ml,    "SDM (SLM ML Durbin)",    wmat)
mr_slm_wx_2sls <- safe_moran_resid(slm_wx_2sls,  "SDM 2SLS",               wmat)
mr_slm_wx_white<- safe_moran_resid(slm_wx_white, "SDM White",              wmat)
mr_sarar_wx_ml <- safe_moran_resid(sarar_wx_ml,  "SARAR ML Durbin",        wmat)


# --- 2. LM TESTS: No WX vs With WX -------------------------------------------
# Compare adjRSerr vs adjRSlag across both specs.
# If adjRSlag becomes significant once WX is included, the earlier SEM
# conclusion may have been partly an artefact of the restricted specification.

cat("\n=== 2. LM TESTS: No WX vs With WX ===\n")
cat("Key: does adjRSlag become significant when WX is added?\n")
if (!is.null(lm_tests))    { cat("\n--- Baseline (no WX) ---\n");   summary(lm_tests) }
if (!is.null(lm_tests_wx)) { cat("\n--- Durbin (with WX) ---\n");   summary(lm_tests_wx) }


# --- 3. AIC / BIC HORSE-RACE (ML models) -------------------------------------
# Lower AIC/BIC = better fit penalised for complexity.
# dAIC > 10 = decisive evidence against that model.

cat("\n=== 3. AIC / BIC MODEL COMPARISON (ML estimators) ===\n")
cat("dAIC > 10 = strong evidence against that model\n\n")

ml_fit <- function(model, label) {
  if (is.null(model)) return(NULL)
  tryCatch(
    data.frame(Model  = label,
               LogLik = round(as.numeric(logLik(model)), 2),
               k      = attr(logLik(model), "df"),
               AIC    = round(AIC(model), 2),
               BIC    = round(BIC(model), 2)),
    error = function(e) { message("SKIPPED ml_fit [", label, "]: ", conditionMessage(e)); NULL }
  )
}

ml_rows <- Filter(Negate(is.null), list(
  ml_fit(ols_lm,      "OLS"),
  ml_fit(ols_wx_lm,   "OLS Durbin (SLX)"),
  ml_fit(sem_ml,      "SEM ML"),
  ml_fit(sem_wx_ml,   "SEM ML Durbin"),
  ml_fit(slm_ml,      "SLM ML"),
  ml_fit(slm_wx_ml,   "SDM (SLM ML Durbin)"),
  ml_fit(sarar_ml,    "SARAR ML"),
  ml_fit(sarar_wx_ml, "SARAR ML Durbin")
))

if (length(ml_rows) > 0) {
  ml_table      <- do.call(rbind, ml_rows)
  ml_table$dAIC <- round(ml_table$AIC - min(ml_table$AIC), 2)
  ml_table$dBIC <- round(ml_table$BIC - min(ml_table$BIC), 2)
  ml_table      <- ml_table[order(ml_table$AIC), ]
  print(ml_table, row.names = FALSE)
}

## Likelihood Ratio Tests (nested comparisons)
cat("\n--- Likelihood Ratio Tests ---\n")
safe_lr(ols_lm,    sem_ml,      "OLS vs SEM ML")
safe_lr(ols_lm,    slm_ml,      "OLS vs SLM ML")
safe_lr(ols_lm,    ols_wx_lm,   "OLS vs OLS Durbin (SLX)")
safe_lr(sem_ml,    sem_wx_ml,   "SEM ML vs SEM ML Durbin")
safe_lr(slm_ml,    slm_wx_ml,   "SLM ML vs SDM")
safe_lr(sem_ml,    sarar_ml,    "SEM ML vs SARAR ML")
safe_lr(slm_ml,    sarar_ml,    "SLM ML vs SARAR ML")
safe_lr(slm_wx_ml, sarar_wx_ml, "SDM vs SARAR ML Durbin")


# --- 4. GMM FIT COMPARISON (RSS, Pseudo-R2, sigma2) --------------------------
# For GMM estimators where AIC is unavailable.
# WARNING: RSS is NOT comparable across models that differ in whether
# they include Wy — SLM will always have lower RSS by construction.

cat("\n=== 4. GMM FIT COMPARISON ===\n")
cat("WARNING: RSS is not comparable across models with and without Wy\n\n")

tss <- sum((data$CHC - mean(data$CHC))^2)
n   <- nrow(data)

gmm_fit2 <- function(model, label) {
  if (is.null(model)) return(NULL)
  tryCatch({
    resid <- as.numeric(residuals(model))
    rss   <- sum(resid^2)
    data.frame(Model     = label,
               RSS       = round(rss, 3),
               Pseudo.R2 = round(1 - rss / tss, 4),
               sigma2    = round(rss / (n - length(model$coef)), 4))
  }, error = function(e) { message("SKIPPED gmm_fit2 [", label, "]: ", conditionMessage(e)); NULL })
}

lm_fit2 <- function(lm_model, label) {
  if (is.null(lm_model)) return(NULL)
  resid <- as.numeric(residuals(lm_model))
  rss   <- sum(resid^2)
  data.frame(Model     = label,
             RSS       = round(rss, 3),
             Pseudo.R2 = round(1 - rss / tss, 4),
             sigma2    = round(rss / (n - length(coef(lm_model))), 4))
}

cat("\n-- Baseline (no WX) --\n")
gmm_rows_base <- Filter(Negate(is.null), list(
  lm_fit2(ols_lm,     "OLS"),
  gmm_fit2(sem_gls,   "SEM KP (FGLS)"),
  gmm_fit2(sem_kp,    "SEM KP robust [preferred]"),
  gmm_fit2(slm_2sls,  "SLM 2SLS"),
  gmm_fit2(slm_white, "SLM White")
))
if (length(gmm_rows_base) > 0) {
  gmm_base <- do.call(rbind, gmm_rows_base)
  gmm_base <- gmm_base[order(gmm_base$RSS), ]
  print(gmm_base, row.names = FALSE)
}

cat("\n-- Durbin (with WX) --\n")
gmm_rows_wx <- Filter(Negate(is.null), list(
  lm_fit2(ols_wx_lm,    "OLS Durbin (SLX)"),
  gmm_fit2(sem_wx_gls,  "SEM KP Durbin (FGLS)"),
  gmm_fit2(sem_wx_kp,   "SEM KP robust Durbin"),
  gmm_fit2(slm_wx_2sls, "SDM 2SLS"),
  gmm_fit2(slm_wx_white,"SDM White")
))
if (length(gmm_rows_wx) > 0) {
  gmm_wx <- do.call(rbind, gmm_rows_wx)
  gmm_wx <- gmm_wx[order(gmm_wx$RSS), ]
  print(gmm_wx, row.names = FALSE)
}


# --- 5. JOINT MODEL SELECTION SUMMARY TABLE -----------------------------------
# All criteria in one place for a level-playing-field comparison.
# A good model: resolves residual SA + competitive AIC + good Pseudo-R2

cat("\n=== 5. JOINT MODEL SELECTION SUMMARY ===\n\n")

moran_p <- function(mt) if (is.null(mt)) NA else round(mt$p.value, 4)

models_summary <- data.frame(
  Model = c(
    "OLS", "OLS Durbin (SLX)",
    "SEM ML", "SEM ML Durbin",
    "SEM KP [preferred]", "SEM KP Durbin",
    "SLM ML", "SDM (SLM ML Durbin)",
    "SLM 2SLS", "SDM 2SLS",
    "SARAR ML", "SARAR ML Durbin"
  ),
  WX = c("No","Yes","No","Yes","No","Yes","No","Yes","No","Yes","No","Yes"),
  Spatial.Param = c(
    "—","—",
    "λ","λ","λ","λ",
    "ρ","ρ","ρ","ρ",
    "ρ+λ","ρ+λ"
  ),
  Resid.Moran.p = c(
    moran_p(mr_ols),       moran_p(mr_ols_wx),
    moran_p(mr_sem_ml),    moran_p(mr_sem_wx_ml),
    moran_p(mr_sem_kp),    moran_p(mr_sem_wx_kp),
    moran_p(mr_slm_ml),    moran_p(mr_slm_wx_ml),
    moran_p(mr_slm_2sls),  moran_p(mr_slm_wx_2sls),
    moran_p(mr_sarar_ml),  moran_p(mr_sarar_wx_ml)
  ),
  stringsAsFactors = FALSE
)

models_summary$SA.Resolved <- ifelse(
  is.na(models_summary$Resid.Moran.p), "N/A",
  ifelse(models_summary$Resid.Moran.p >= 0.05, "YES", "NO *")
)

print(models_summary, row.names = FALSE)
cat("\n* NO  = residual spatial autocorrelation remains (p < 0.05) — model inadequate\n")
cat("  YES = residuals spatially random — spatial structure adequately captured\n")
cat("\nConclusion: the preferred model is the one that achieves SA.Resolved = YES\n")
cat("with the lowest AIC (among ML models) or best Pseudo-R2 (among GMM models),\n")
cat("and has a theoretically defensible story for the rise of Chega.\n")


# =============================================================================
# DIRECT / INDIRECT / TOTAL EFFECTS
# =============================================================================

extract_effects <- function(model, model_name = "") {
  if (is.null(model)) {
    message("SKIPPED effects [", model_name, "]: model not available")
    return(invisible(NULL))
  }
  tryCatch({
    if (inherits(model, "Sarlm")) {
      rho     <- if (!is.null(model$rho)) model$rho else 0
      beta    <- model$coefficients[-1]
      beta_se <- model$rest.se[-1]
    } else if (inherits(model, "sphet")) {
      coefs  <- model$coef
      se_all <- sqrt(diag(model$var))
      mtype  <- tolower(as.character(model$model)[1])
      if (mtype == "ols") {
        rho <- 0; beta <- coefs[-1]; beta_se <- se_all[-1]
      } else if (mtype %in% c("lag","ivhac")) {
        rho <- tail(coefs,1); beta <- head(coefs,-1)[-1]; beta_se <- head(se_all,-1)[-1]
      } else if (mtype == "error") {
        rho <- 0; beta <- head(coefs,-1)[-1]; beta_se <- head(se_all,-1)[-1]
      } else if (mtype == "sarar") {
        rho <- tail(coefs,2)[1]; beta <- head(coefs,-2)[-1]; beta_se <- head(se_all,-2)[-1]
      } else {
        warning(paste("extract_effects: unrecognised model family:", mtype))
        return(invisible(NULL))
      }
    } else {
      warning(paste("extract_effects: unrecognised class for", model_name))
      return(invisible(NULL))
    }
    eff <- data.frame(
      Direct   = round(beta,                              4),
      Indirect = round((rho/(1-rho)) * beta,              4),
      Total    = round(beta + (rho/(1-rho)) * beta,       4),
      Std.Err  = round(beta_se,                           4),
      p.value  = round(2*(1-pnorm(abs(beta/beta_se))),    4),
      row.names = names(beta)
    )
    cat("\n===", model_name, "===\n")
    print(eff)
    return(invisible(eff))
  }, error = function(e) {
    message("SKIPPED effects [", model_name, "]: ", conditionMessage(e))
    invisible(NULL)
  })
}

## OLS
extract_effects(ols,          "OLS")
extract_effects(ols_wx,       "OLS Durbin (SLX)")
extract_effects(ols_white,    "OLS White")
extract_effects(ols_wx_white, "OLS White Durbin")

## SEM ML
extract_effects(sem_ml,    "SEM ML")
extract_effects(sem_wx_ml, "SEM ML Durbin")

## SEM KP
extract_effects(sem_gls,    "SEM KP (FGLS)")
extract_effects(sem_wx_gls, "SEM KP Durbin (FGLS)")
extract_effects(sem_kp,     "SEM KP robust [preferred]")
extract_effects(sem_wx_kp,  "SEM KP robust Durbin")

## SLM
extract_effects(slm_ml,       "SLM ML")
extract_effects(slm_wx_ml,    "SDM (SLM ML Durbin)")
extract_effects(slm_2sls,     "SLM 2SLS")
extract_effects(slm_wx_2sls,  "SDM 2SLS")
extract_effects(slm_white,    "SLM White")
extract_effects(slm_wx_white, "SDM White")

## SARAR
extract_effects(sarar_ml,    "SARAR ML")
extract_effects(sarar_wx_ml, "SARAR ML Durbin")

# =============================================================================
# END OF SCRIPT
# =============================================================================
