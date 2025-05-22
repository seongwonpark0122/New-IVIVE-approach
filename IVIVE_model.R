### Load library
library(dplyr)
library(nlmixr2)

### Load observed data
df = …

### Save model as function
model = function() {
  ini({
    ### Set the initial value for each parameters
    ### Weibull distribution parameters for parent drug and metabolite
    lvmax_p <- log(0.0631)
    lvmax_m <- log(0.03)
    alpha <- fix(1.565)
    beta <- fix(-3.3)
    
    lmeta <- log(4.71e-06) # cellular metabolic clearance (mL/min/cell)
    lportion <- log(0.0155) # metabolic conversion fraction
    
    V2 <- fix(1) # volume of insert = 1 mL
    
    sf <- 0.643 # scaling factor for cell counts
    
    eta.meta ~ fix(0) # fixed variability as 0 in model development step
    # eta.meta can be changed in simulation steps
    
    add.err <-  fix(0) # 0% for random variability of parent drug
    add.err_M <-  fix(0) # 0% for random variability of metabolite
  })
  model({
    ### Change model parameters from log scale to linear scale
    vmax_p = exp(lvmax_p)
    vmax_m = exp(lvmax_m)
    portion = expit(lportion)    
    
    ### Designate rate constants for parent drug and metabolite with Weibull distribution equation
    wb = (1 - exp(-((time / alpha)^beta)))
    
    ### Determine metabolic clearance by incorporating the influence of cell count.
    ### Cell count information is included in the input data. 
    ### Alternatively, this information can be directly entered in the initialization section.
    CL_meta = exp(lmeta + sf*log(Cell_Count) + eta.meta)
    elimA2 = CL_meta * A2    
    
    ### Save flow rates for each compartment
    flowA = vmax_p * (A1 - A2) * wb # parent flow
    flowAm = vmax_m * (A1_Met - A2_Met) * wb # metabolite flow    
    
    ### Differential equation for each compartments
    d/dt(A1) = - flowA ### media compartment of parent drug
    d/dt(A2) = flowA - elimA2 ### insert compartment of parent drug
    d/dt(A1_Met) = - flowAm ### media compartment of parent drug
    d/dt(A2_Met) = flowAm + elimA2 * portion * (312.1/296.148) ### insert compartment of metabolite
    
    ### Calculate concentration for parent drug and metabolite
    ### Specifies concentration for observed values in linear scale, or log-transformed concentration for observed values in log scale.
    ### C2 is for concentration of parent drug in insert, and C2_Met is for the metabolite
    C2 <- log(A2/V2)
    C2_Met <- log(A2_Met/V2)
    
    C2 ~ add(add.err) | A2
    C2_Met ~ add(add.err_M) | A2_Met
  })
}

### Comfile the model
mod <- nlmixr2(model)

### Conduct estimation
fit <- nlmixr2(mod, df, "focei", control = foceiControl(seed = 1234))
