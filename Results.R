# Colorado Emission Results (2015-2020)
CO_results <- ElectricityconsumptionEmissions(input_filename = "Colorado Input Data.csv", 
                                              emission_source = 12, nyears = 6, 
                                              iseed = 131313, nreps = 10000)

# Vermont Emission Results (2015-2020)
VT_results <- ElectricityconsumptionEmissions(input_filename = "Vermont Input Data.csv", 
                                              emission_source = 12, nyears = 6, 
                                              iseed = 131313, nreps = 10000)