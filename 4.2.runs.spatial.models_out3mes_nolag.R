source("4.1.load.packages.data_out3mes_nolag.R")



baseformula <- Y ~ 1 + basis_RR + basis_Precip + basis_Nino34SSTA + basis_TNA + basis_NDVI + basis_LSD +
  f(Month, replicate = nCanton, model = "rw1", cyclic = TRUE, constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior)

formula.iid <- update.formula(baseformula, ~. + f(nCanton, model="iid"))

formula.1.1 <- update.formula(baseformula, ~. + f(nCanton, model = "besag", graph = W.canton32))
formula.1.2 <- update.formula(baseformula, ~. + f(nCanton, model = "besagproper", graph = W.canton32))
formula.1.3 <- update.formula(baseformula, ~. + f(nCanton, model = "bym", graph = W.canton32))
formula.1.4 <- update.formula(baseformula, ~. + f(nCanton, model = "bym2", graph = W.canton32))

formula.2.1 <- update.formula(baseformula, ~. + f(nCanton, model = "besag", graph = W.adj1))
formula.2.2 <- update.formula(baseformula, ~. + f(nCanton, model = "besagproper", graph = W.adj1))
formula.2.3 <- update.formula(baseformula, ~. + f(nCanton, model = "bym", graph = W.adj1))
formula.2.4 <- update.formula(baseformula, ~. + f(nCanton, model = "bym2", graph = W.adj1))

formula.3.1 <- update.formula(baseformula, ~. + f(nCanton, model = "besag", graph = W.adj2))
formula.3.2 <- update.formula(baseformula, ~. + f(nCanton, model = "besagproper", graph = W.adj2))
formula.3.3 <- update.formula(baseformula, ~. + f(nCanton, model = "bym", graph = W.adj2))
formula.3.4 <- update.formula(baseformula, ~. + f(nCanton, model = "bym2", graph = W.adj2))

formula.1.1a <- update.formula(baseformula, ~. + f(nCanton, model = "besag",  replicate = nYear, graph = W.canton32))
formula.1.2a <- update.formula(baseformula, ~. + f(nCanton, model = "besagproper",  replicate = nYear, graph = W.canton32))
formula.1.3a <- update.formula(baseformula, ~. + f(nCanton, model = "bym",  replicate = nYear, graph = W.canton32))
formula.1.4a <- update.formula(baseformula, ~. + f(nCanton, model = "bym2",  replicate = nYear, graph = W.canton32))

formula.2.1a <- update.formula(baseformula, ~. + f(nCanton, model = "besag",  replicate = nYear, graph = W.adj1))
formula.2.2a <- update.formula(baseformula, ~. + f(nCanton, model = "besagproper",  replicate = nYear, graph = W.adj1))
formula.2.3a <- update.formula(baseformula, ~. + f(nCanton, model = "bym",  replicate = nYear, graph = W.adj1))
formula.2.4a <- update.formula(baseformula, ~. + f(nCanton, model = "bym2",  replicate = nYear, graph = W.adj1))

formula.3.1a <- update.formula(baseformula, ~. + f(nCanton, model = "besag",  replicate = nYear, graph = W.adj2))
formula.3.2a <- update.formula(baseformula, ~. + f(nCanton, model = "besagproper",  replicate = nYear, graph = W.adj2))
formula.3.3a <- update.formula(baseformula, ~. + f(nCanton, model = "bym",  replicate = nYear, graph = W.adj2))
formula.3.4a <- update.formula(baseformula, ~. + f(nCanton, model = "bym2",  replicate = nYear, graph = W.adj2))


mymodel <- function(formula, data = df1, family = "nbinomial",OFF= OFF , config = FALSE)
{
  model <- inla(formula = formula, data = data, family = family, offset = OFF,
                control.inla = list(strategy = 'adaptive'), 
                control.compute = list(dic = TRUE, config = config, 
                                       cpo = TRUE, return.marginals = FALSE),
                control.fixed = list(correlation.matrix = TRUE, 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE), 
                verbose = FALSE)
  model <- inla.rerun(model)
  return(model)
}

# create a list of formulas
formulas <- list(formula.iid, formula.1.1, formula.1.2, formula.1.3, formula.1.4,
                 formula.2.1, formula.2.2, formula.2.3, formula.2.4,
                 formula.3.1, formula.3.2, formula.3.3, formula.3.4, 
                 formula.1.1a, formula.1.2a, formula.1.3a, formula.1.4a,
                 formula.2.1a, formula.2.2a, formula.2.3a, formula.2.4a,
                 formula.3.1a, formula.3.2a, formula.3.3a, formula.3.4a)
# create model label string
lab <- c("formula.iid", "formula.1.1", "formula.1.2", "formula.1.3", "formula.1.4",
         "formula.2.1", "formula.2.2", "formula.2.3", "formula.2.4",
         "formula.3.1", "formula.3.2", "formula.3.3", "formula.3.4", 
         "formula.1.1a", "formula.1.2a", "formula.1.3a", "formula.1.4a",
         "formula.2.1a", "formula.2.2a", "formula.2.3a", "formula.2.4a",
         "formula.3.1a", "formula.3.2a", "formula.3.3a", "formula.3.4a")

# create a function to run a model for each formula in the list and save the model output to file
# WARNING: this may take a long time to run
models <- lapply(1:length(formulas), 
                 function(i) {
                   print(i)
                   model <- mymodel(formulas[[i]], df1)
                   save(model, file = paste0("output_out3mes_nolags/", lab[i],".RData"))})






