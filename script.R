#R script for calculating the mass balance of in-situ marine methane concentration.
#D'Angelo A. and Loose B. (2020) - University of Rhode Island, Graduate School of Oceanography

#load dataframe "data"

library(dplyr)

#Bunsen solubility coefficient: https://github.com/URIGSO/Gas-Solubility-Codes/tree/R_code

#calculate the methane gas solubility according to the Bunsen solubility coefficient (Yamamoto et al., 1976)
    bunsen <- function(T, S){
    
    A1 <- -67.1962
    A2 <- 99.1624
    A3 <- 27.9015
    B1 <- -0.072909
    B2 <- 0.041674
    B3 <- -0.0064603
    T<- 1
    T.K <- T + 273.15
    
    bunsen.coeff.L.L.x <- exp(A1 + A2*(100/T.K) + A3*(log(T.K/100)) + S*(B1 + B2*(T.K/100) + B3*((T.K/100)^2)))
    
    bunsen.coeff.L.L <- c(bunsen.coeff.L.L.x)
    
    return(bunsen.coeff.L.L)
  }
  #Calculate Bunsen solubility coefficients using water bath T and in-situ S
  beta_sol=bunsen(10,data$Sal)
  
  # Ideal gas constant (mol/L)
  IGC=(101.325/(8.314*263))  #101.325/(8.314*(273-Celsius)) 
  print(IGC)
  
  # Pressure correction
  P_SSIM=data$Pssim/(data$airP/1.3332) #unitless
  
  # Volume correction
  a=1027 #g/l as eight of bags is in g
  Vw=MOSAiC.Ambient.calibrated$bag.M/a #L
  print(Vw)
  
  # Dilution factor correction
  data$vol.L<- data$vol.mL/1000 #L
  SSIM.dil.factor=(0.022/data$vol.L) 
  
  # Calculate methane moles in headspace
  #pCH4 is the partial pressure of methane recorded from the Picarro analyzer
  pCH4.hs=SSIM.dil.factor*P_SSIM*data$pCH4 #ppm
  chi.hs=pCH4.hs*10^-6 #ppp
  beta_sol.moles=beta_sol*IGC #mol/L
  print(beta_sol.moles)
  
  # Calculate the methane concentration in headspace
  cCH4.hs=beta_sol.moles*chi.hs #mol/L
  
  # Calculate the mass of gas in the water
  Mg_w=Vw*cCH4.hs #mol
  
  # Calculate the mass of gas in the headspace
  #Vg=Vhs in ambient bags
  Vg=0.05 #l
  Mg_hs=Vg*chi.hs*IGC  #mol
  
  # Calculate the total methane moles
  data$Mtot=Mg_w+Mg_hs #mol
  
  # Calculate methane concentration in the samples (nmol/L)
  data$CH4_conc <- (data$Mtot/Vw)*10^9 
  
  # Calculate methane saturation with respect to atmospheric methane  
  bunsen <- function(T, S){
    
    A1 <- -67.1962
    A2 <- 99.1624
    A3 <- 27.9015
    B1 <- -0.072909
    B2 <- 0.041674
    B3 <- -0.0064603
    T<- 1
    T.K <- T + 273.15
    
    bunsen.coeff.L.L.x <- exp(A1 + A2*(100/T.K) + A3*(log(T.K/100)) + S*(B1 + B2*(T.K/100) + B3*((T.K/100)^2)))
    
    bunsen.coeff.L.L <- c(bunsen.coeff.L.L.x)
    
    return(bunsen.coeff.L.L)
  }
  
  beta_sol_eq=bunsen(data$potT,data$Sal) #in-situ T and S
  beta_sol_eq_moles = beta_sol_eq*IGC # mol/L
  
  # Set amtospheric values for methane
  Chi_atm = 1.8/10^6  #ppp
  # Calculate equilibrium concentration
  data$Conc_eq = (beta_sol_eq_moles*Chi_atm)*10^9 #nmol/L
  print(data$Conc_eq)
  
  # Calculate percentages
  data$Conc_eq_perc =(data$CH4_conc/data$Conc_eq)*100 #%
  
  # Calculate saturation anomaly: (C_insitu/C_eq -1)*100
  data$sat.anomaly<- (data$CH4_conc/data$Conc_eq -1)*100
  print(data$sat.anomaly)
  
  #Reference: Yamamoto, S., et al. 1976. Solubility of methane in distilled water and seawater. Journal of Chemical & Engineering Data, 21, 78-80.
  
