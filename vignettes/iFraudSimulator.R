## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ---- echo = FALSE------------------------------------------------------------
suppressMessages(library(iFraudSimulator))
suppressMessages(library(magrittr))
library(bookdown)

## -----------------------------------------------------------------------------
?sfnGenerator

## ---- echo = FALSE, eval = FALSE----------------------------------------------
#  SimObj = sfnGenerator(TargetPrev = 0.05, NrPH = 10000, Seed = 1, printProgress = FALSE, Parallel = FALSE)

## -----------------------------------------------------------------------------
data("SimObj")

## ---- fig.align = 'center', fig.cap = "Summary plot of the synthetic data set", fig.topcaption = TRUE, out.width="100%"----
SimObj

## -----------------------------------------------------------------------------
names(SimObj)

## -----------------------------------------------------------------------------
Dt = copy(SimObj$Dt)

## -----------------------------------------------------------------------------
Homophily(SimObj)

## -----------------------------------------------------------------------------
formulaExpert   = update.formula(SimObj$call$Formulas$Fraud, ExpertJudgement ~ . - Criminal)
formulaExpert   = do.call("substitute", list(expr = formulaExpert, env = alist(AgePH = AgePHScaled)))
netwFeatures = c("fraudScore", colnames(Dt)) %>% .[grepl("fraudScore|n1|n2", .)] %>% .[!grepl("n1Size|n2Size|Orig", .)]
formulaM1 = update.formula(formulaExpert, paste0(". ~ . -", 
                                                 paste0(c(netwFeatures, "n1Size", "n2Size"), collapse = " - ")))
formulaM2 = formulaExpert

## -----------------------------------------------------------------------------
OptFun = list(
  AUC = function(y, p) {
    PRROC::roc.curve(scores.class0 = p, weights.class0 = y)$auc
  },
  TopDecileLift = function(y, p) {
    lift::TopDecileLift(p, y)
  }
)
compModels <- function(SimObj) {
  Dt     = SimObj$Dt
  m1 = glm(formulaM1, family = binomial, data = Dt[Investigated == 1])
  m2 = glm(formulaM2, family = binomial, data = Dt[Investigated == 1])
  
  y    = as.numeric(Dt[Investigated == 0, get("Fraud")])
  p.m1 = predict(m1, newdata = Dt[Investigated == 0], type = "response")
  p.m2 = predict(m2, newdata = Dt[Investigated == 0], type = "response")
  
  Res = 
    list(
      Model1 = list(
        coef = coef(m1), 
        IS   = do.call("c", lapply(OptFun, function(f) f(y = m1$y, p = fitted(m1)))), 
        OOS  = do.call("c", lapply(OptFun, function(f) f(y = y, p = p.m1)))
        ),
      Model2 = list(
        coef = coef(m2), 
        IS   = do.call("c", lapply(OptFun, function(f) f(y = m2$y, p = fitted(m2)))), 
        OOS  = do.call("c", lapply(OptFun, function(f) f(y = y, p = p.m2)))
      )
    )
  return(Res)
}

## -----------------------------------------------------------------------------
compModels(SimObj)

