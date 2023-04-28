"ElectricityconsumptionEmissions" <- 
  function(input_filename = "NAME", emission_source = 1, nyears = 1, iseed = 131313, 
           nreps = 10000)
    # Script developed by K. Peterson
    # Originally developed: April 06, 2023
    # Last upate: April 06, 2023
    #
    ####### Begin Script
  {
    ####### Set seed
    set.seed(iseed)
    #
    ####### Load Library
    library(triangle)
    ####### Import files
    input_data <- read.csv(file=input_filename, header=TRUE, sep=",", fill=FALSE)
    # Ensure all input data for fuel amounts and emission factors are numeric
    for (n in (1:(9+(nyears*2)))){
      input_data[,n+1] <- as.numeric(input_data[,n+1])
    }
    #
    ####### Check Validity of Input Data
    check_fuel_type <- length(input_data[,1])==emission_source
    if(!check_fuel_type){stop("Warning: The number of fuel types provided in the 
                              function call is not equal to the number in the 
                              input file.")}
    #
    # Emission Factors
    for(efact in (1:emission_source)){
      check_CO2_EF <- input_data[efact, 2]>=0
      if (!check_CO2_EF) {stop("Error: Default CO2 emission factor must be greater 
                          than or equal to 0 - Check input file.")}
    }
    check_CO2_EF_min <- input_data[efact, 3] <= input_data[efact,2] & 
      input_data[efact,3] >= 0
    if(!check_CO2_EF_min) {stop("Error: Minimum CO2 emission factor must be 
                          greater than or equal to 0 and less than the default 
                          emission factor vlaue - Check input file.")
    }
    check_CO2_EF_max <- input_data[efact, 4] >= input_data[efact,2]
    if(!check_CO2_EF_max) {stop("Error: Maximum CO2 emission factor must be 
                          greater than or equal to the default emission factor
                          value - Check input file.")
    }
    check_CH4_EF <- input_data[efact, 5] >= 0
    if (!check_CH4_EF) {stop("Error: Default CH4 emission factor must be greater 
                        than or equal to 0 - Check input file.")
    }
    check_CH4_EF_min <- input_data[efact, 6] <= input_data[efact,5] & 
      input_data[efact,6] >= 0
    if(!check_CH4_EF_min) {stop("Error: Minimum CH4 emission factor must be 
                          greater than or equal to 0 and less than the default 
                          emission factor vlaue - Check input file.")
    }
    check_CH4_EF_max <- input_data[efact, 7] >= input_data[efact,5]
    if(!check_CH4_EF_max) {stop("Error: Maximum CH4 emission factor must be 
                          greater than or equal to the default emission factor 
                          value - Check input file.")
    }
    check_N2O_EF <- input_data[efact, 8] >= 0
    if (!check_N2O_EF) {stop("Error: Default N2O emission factor must be greater 
                        than or equal to 0 - Check input file.")
    }
    check_N2O_EF_min <- input_data[efact, 9] <= input_data[efact,8] & 
      input_data[efact,9] >= 0
    if(!check_N2O_EF_min) {stop("Error: Minimum N2O emission factor must be greater 
                          than or equal to 0 and less than the default emission 
                          factor vlaue - Check input file.")
    }
    check_N2O_EF_max <- input_data[efact, 10] >= input_data[efact,8]
    if(!check_N2O_EF_max) {stop("Error: Maximum N2O emission factor must be 
                          greater than or equal to the default emission factor 
                          value - Check input file.")
    }
    #
    # Usage Amounts
    for(y in (1:nyears)) {
      for (f in (1:emission_source)){
        check_usage <- input_data[f,11+((y-1)*2)] >= 0
        if(!check_usage) {stop("Error: Usage must be greater than or equal to 0 
                          - Check input file.")
        }
        check_usage_sd<-input_data[f,12+((y-1)*2)] >= 0
        if(!check_usage_sd) {stop("Error: Usage standard deviation must be greater 
                            than or equal to 0 - Check input file.")
        }
      }
    }
    #
    ####### Monte Carlo Analysis 
    ##### Simulate nreps of emission factors (triangle distribution) 
    #     and usage amounts (normal distribution)
    ### Emission Factors
    # CO2
    CO2_EF_sim <- matrix(0,nrow=emission_source, ncol=nreps)
    for (f in (1:emission_source)){
      CO2_EF_sim[f,]<-rtriangle(nreps, a=input_data[f,3],b=input_data[f,4],
                                c=input_data[f,2])
    }
    #
    # CH4
    CH4_EF_sim <- matrix(0,nrow=emission_source, ncol=nreps)
    for (f in (1:emission_source)){
      CH4_EF_sim[f,]<-rtriangle(nreps, a=input_data[f,6],b=input_data[f,7],
                                c=input_data[f,5])
    }
    #
    # N2O
    N2O_EF_sim <- matrix(0,nrow=emission_source, ncol=nreps)
    for (f in (1:emission_source)){
      N2O_EF_sim[f,]<-rtriangle(nreps, a=input_data[f,9],b=input_data[f,10],
                                c=input_data[f,8])
    }
    #
    ### Usage Amounts
    usage_sim<-matrix(0,nrow=emission_source*nyears,ncol=nreps)
    for (y in (1:nyears)){
      for(f in(1:emission_source)){
        usage_sim[f+(emission_source*(y-1)),]<-rnorm(nreps,mean=input_data[f,11+((y-1)*2)],
                                                     sd=input_data[f,12+((y-1)*2)])
      }
    }
    #
    ####### Calculate Emissions 
    # IPCC 2006 GL: Emissions = Usage * EF
    # Units: Emissions in MT, usage in MWh, and emission factors in MT/MWh
    ### Deterministic Estimation
    Deterministic_CO2eq<-matrix(0,nrow=emission_source,ncol=nyears)
    for(y in (1:nyears)){
      for (d in (1:emission_source)){
        Deterministic_CO2eq[d,y]<-(input_data[d,2]*input_data[d,11+((y-1)*2)])+
          (input_data[d,5]*input_data[d,11+((y-1)*2)])*28+
          (input_data[d,8]*input_data[d,11+((y-1)*2)])*265
      }
    }
    #
    ## Total CO2eq emissions for each year
    Deterministic_CO2eq_total <- apply(Deterministic_CO2eq, MAR=2, FUN="sum")
    #  
    ### Probabilistic Estimation
    Probabilistic_CO2eq <- matrix(0, nrow=emission_source*nyears,ncol=nreps)
    for (y in (1:nyears)){
      for(p in (1:emission_source)){
        Probabilistic_CO2eq[p+(emission_source*(y-1)),]<-
          ((CO2_EF_sim[p,]*usage_sim[p+(emission_source*(y-1)),]))+
          ((CH4_EF_sim[p,]*usage_sim[p+(emission_source*(y-1)),])*28)+
          ((N2O_EF_sim[p,]*usage_sim[p+(emission_source*(y-1)),])*265)
      }
    }
    #
    ## Total CO2eq emissions
    Probabilistic_CO2eq_total <- matrix(0,nrow=nyears, ncol=nreps)
    for (y in(1:nyears)){
      Probabilistic_CO2eq_total[y,]<-apply(Probabilistic_CO2eq[(1+((y-1)*emission_source)):
                                                                 (emission_source+((y-1)*emission_source)),],
                                           MAR=2, FUN="sum")
    }
    #
    ####### Estimate median and confidence intervals
    # Create matrix with Col 1 = median, col 2 = 2.5 percentile, 
    # and col 3 = 97.5 percentile
    emission_results <- matrix (0, nrow=nyears, ncol=3)
    for (y in (1:nyears)) {
      emission_results[y,1] <- median(Probabilistic_CO2eq_total[y,])
      q <- quantile(Probabilistic_CO2eq_total[y,], probs = c(0.025, 0.975))
      emission_results[y, 2] <- q[1]
      emission_results[y, 3] <- q[2]
      check_emissions <- Deterministic_CO2eq_total[y] >= emission_results[y,2] &
        Deterministic_CO2eq_total[y] <= emission_results[y,3]
      if(!check_emissions) {
        cat("WARNING: Deterministic Solution for year", y, "is outside of its respective 
            confidence interval.")
      }
    }
    #
    ####### Return Statement 
    # Convert to MMT CO2 from MT CO2
    emission_results <- emission_results/10^6
    colnames(emission_results) <- c("Median MMTCO2eq", "2.5 Percentile", "97.5 Percentile")
    return (emission_results)
    #
    ####### End Script
  }
