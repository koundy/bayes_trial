
library(MCMCglmm)
library(lme4)
library(ggplot2)
library(reshape2)
#library(brms)
source("~/ggplot_theme_Publication-2.R")
theme_set(theme_Publication())

#########################  Generation of random dataset 

set.seed(102)
n_indiv <- 200
n_per_indiv <- 4
n_tot <- n_indiv * n_per_indiv

id <- gl(n_indiv, n_per_indiv) ## generate identifier columns

av_wealth <- rlnorm(n_indiv, 0, 1) ## avg wealth of individuals
ac_wealth <- av_wealth[id] + rlnorm(n_tot, 0, 1) ## wealth in each year:

av_ratio <- rbeta(n_indiv, 10, 10) ## 50/50, close to Gaussian
## car/holiday ratio by year (larger shape parameter/less variability)
ac_ratio <- rbeta(n_tot, 2 * av_ratio[id], 2 * (1 - av_ratio[id]))

y.car <- (ac_wealth * ac_ratio)^0.25
y.hol <- (ac_wealth * (1 - ac_ratio))^0.25
Spending <- data.frame(y.hol, y.car, id)

###########################  Visualization of data 

qplot(y.hol, y.car, data = Spending) + geom_smooth(method = "lm")

m0 <- lm(y.car ~ y.hol, data = Spending)

m1_id <- MCMCglmm(y.car ~ y.hol, random = ~id, data = Spending, verbose = FALSE)
m1_idyr <- MCMCglmm(cbind(y.hol, y.car) ~ trait - 1,
                    random = ~ us(trait):id,
                    rcov = ~ us(trait):units, data = Spending, family = c("gaussian", "gaussian"),
                    verbose = FALSE
)

rd <- as.data.frame(m1_idyr$VCV)
id.regression <- with(rd, `y.car:y.hol.id` / `y.hol:y.hol.id`)
units.regression <- with(rd, `y.car:y.hol.units` / `y.hol:y.hol.units`)
res_m1_id <- setNames(summary(m1_id)$solutions[2, 1:3], c("est", "lwr", "upr"))
res_m1_idyr_id <- c(est = mean(id.regression), setNames(quantile(
    id.regression,
    c(0.025, 0.975)
), c("lwr", "upr")))
res_m1_idyr_units <- c(est = mean(units.regression), setNames(quantile(
    units.regression,
    c(0.025, 0.975)
), c("lwr", "upr")))

ff <- function(x) factor(x, levels = unique(x))
combdat <- data.frame(
    param = ff(c("naive", "m_L1", "m_L2.units", "m_L2.id")),
    type = ff(c("lm", rep("MCMCglmm", 3))), rbind(
        res_reg0, res_m1_id, res_m1_idyr_units,
        res_m1_idyr_id
    )
)


(g1 <- ggplot(combdat, aes(param, est, ymin = lwr, ymax = upr, colour = type)) +
        geom_pointrange() +
        labs(x = "", y = "estimate") +
        geom_hline(
            yintercept = 0,
            colour = "black", lwd = 2, alpha = 0.2
        ) +
        coord_flip() +
        scale_colour_brewer(palette = "Dark2"))


mSpending <- melt(data.frame(Spending, obs = seq(nrow(Spending))), id.var = c(
    "obs",
    "id"
), variable.name = "trait")

L1_id <- lmer(y.car ~ y.hol + (1 | id), data = Spending)

L2_unid <- lmer(value ~ trait - 1 + (0 + trait | id) + (0 + trait | obs), data = mSpending)

vc_unid <- VarCorr(L2_unid)
## transform VarCorr output into regression coefficients
vchack <- function(vc, resvar = NULL) {
    vc1 <- vc$obs
    if (is.null(resvar)) {
        resvar <- vc1[1, 1]
    }
    diag(vc1) <- diag(vc1) + resvar
    res_idyr <- setNames(c(vc1["traity.hol", "traity.car"] / vc1[
        "traity.hol",
        "traity.hol"
    ], NA, NA), c("est", "lwr", "upr"))
    vc2 <- vc$id
    res_id <- setNames(c(
        vc2["traity.hol", "traity.car"] / vc2["traity.hol", "traity.hol"],
        NA, NA
    ), c("est", "lwr", "upr"))
    rbind(idyr = res_idyr, id = res_id)
}
vchack(vc_unid, resvar = sigma(L2_unid)^2)

f0 <- lmer(value ~ trait - 1 + (0 + trait | id) + (0 + trait | obs), data = mSpending)
names(getME(f0, "theta"))

tmpf <- lmer(value ~ trait - 1 + (0 + trait | id) + (0 + trait | obs),
             data = mSpending,
             devFunOnly = TRUE
)
tmpf2 <- function(theta2) {
    tmpf(c(1, theta2))
}


library(optimx)
opt1 <- optimx(par = c(0, 1, 1, 0, 1), fn = tmpf2, method = "bobyqa")

theta <- c(1, opt1$par$par) ## estimated Cholesky factors
pwrss <- opt1$fvalues$fvalues ## penalized weighted residual sum of sq
n <- nrow(mSpending)
cnms <- f0@cnms ## named list of names of variables per RE
vc <- lme4:::mkVarCorr(
    sc = pwrss / n,
    cnms = cnms,
    nc = sapply(cnms, length), ## vars per random effect term
    theta = theta,
    nms = names(cnms)
)
attr(vc, "useSc") <- TRUE
class(vc) <- "VarCorr.merMod"
vc

lme4res <- vchack(vc)
combdat2 <- rbind(combdat, data.frame(
    param = c("l_L2.units", "l_L2.id"), type = "lme4",
    vchack(vc)
))

g1 %+% combdat2
