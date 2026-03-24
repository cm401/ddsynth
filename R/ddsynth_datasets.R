
datasets_Nipah <- list(
  d1  = list(median = 10.0, min =  9.0, max = 12.0, n =  4, source = ""),
  d2  = list(median =  4.0, min =  2.0, max =  7.0, n =  6, source = ""),
  d3  = list(median =  9.0, min =  6.0, max = 14.0, n = 11, source = ""),
  d4  = list(median =  9.5, min =  4.0, max = 14.0, n = 22, source = ""),
  d5  = list(median =  8.0, min =  3.0, max = 20.0, n = 15, source = ""),
  d6  = list(median =  9.0, min =  6.0, max = 11.0, n = 14, source = ""),
  d7  = list(median = 10.0, min =  8.0, max = 15.0, n = 11, source = ""),
  d8  = list(median =  9.0, min =  6.0, max = 11.0, n = 11, source = ""),
  d9  = list(median =  9.0, Q1  =  8.0, Q3  = 11.0, n = 82, source = ""),
  d10 = list(mean   =  9.3, sd  =  1.9,             n = 18, source = ""),
  d11 = list(
    # n defaults to sum(freq_count) = 11
    freq_value = c(6, 7, 8, 9, 11, 12, 14),
    freq_count = c(1, 1, 3, 2,  1,  2,  1),
    source     = ""
  )
)


# Underlying databased on PERG MVD Systematic Review, published in the Lancet ID
# https://doi.org/10.1016/S1473-3099(23)00515-7

datasets_MVD <- list(         # Marburg Virus Disease
  d1 = list(median =  7.0, min =  2.0, max = 14.0, n = 66, source = ""),
  d2 = list(median = 10.0, Q1  =  8.0, Q3  = 13.0, n = 76, source = "")
)

# Underlying databased on PERG EVD Systematic Review, published in the Lancet ID
# https://doi.org/10.1016/S1473-3099(24)00374-8

datasets_EVD <- list(         # Ebola Virus Disease
  d1  = list(median =  7.0, min =  2.0, max = 20.0, n =  116, source = ""),
  d2  = list(median =  6.0, min =  1.0, max = 16.0, n =   24, source = ""),
  d3  = list(mean   =  9.2, sd  =  6.7,             n =   33, source = ""),
  d4  = list(mean   =  8.6, sd  =  6.1,             n =   20, source = ""),
  d5  = list(mean   =  9.5, sd  =  4.0,             n =   76, source = ""),
  d6  = list(mean   = 10.0, sd  =  1.0,             n =  291, source = ""),  # NOTE: unusually low SD — flag for pre_inference_checks()
  d7  = list(mean   =  9.9, sd  =  5.5,             n =  152, source = ""),
  d8  = list(mean   =  9.7, sd  =  3.7,             n =    8, source = ""),
  d9  = list(mean   =  9.3, sd  =  1.9,             n =   20, source = ""),
  d10 = list(        # West Africa 2014
    # n defaults to sum(freq_count) = 143
    freq_value = c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                   21, 22, 23, 25, 27, 30, 32, 35, 38, 42),
    freq_count = c( 1,  2,  4,  6,  8,  9, 10, 11, 12, 13,
                   11, 10,  8,  7,  5,  4,  3,  3,  2,  2,
                    2,  1,  1,  1,  1,  1,  1,  1,  1,  1),
    source     = ""
  ),
  d11 = list(        # Guinea
    # n defaults to sum(freq_count) = 58
    freq_value = c( 2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
                   12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 30, 38),
    freq_count = c( 1,  2,  3,  4,  4,  4,  4,  5,  5,  4,
                    4,  3,  3,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1),
    source     = ""
  ),
  d12 = list(        # Liberia
    # n defaults to sum(freq_count) = 52
    freq_value = c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 25, 27, 32, 35, 42),
    freq_count = c( 1,  1,  1,  2,  2,  3,  4,  4,  4,  5,
                    4,  4,  3,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1),
    source     = ""
  ),
  d13 = list(        # Sierra Leone
    # n defaults to sum(freq_count) = 30
    freq_value = c( 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
    freq_count = c( 1,  1,  2,  2,  2,  3,  3,  3,  3,  2,  2,  2,  1,  1,  1,  1),
    source     = ""
  )
)

# Underlying databased on PERG Lassa Fever Systematic Review, published in the Lancet Global Health
# https://doi.org/10.1016/S2214-109X(24)00379-6

datasets_Lassa <- list(           # Lassa Fever
  # Single dataset — tau not identifiable; predictions computed at mean(loc_d).
  # See prepare_stan_data_from_datasets() documentation.
  d1 = list(
    # n defaults to sum(freq_count) = 15
    freq_lower = c( 1,  2,  2,  3,  4,  6,  7,  8, 10, 11),
    freq_upper = c(12, 13, 17, 18, 19, 13, 22, 12, 10, 15),
    freq_count = c( 1,  1,  1,  2,  2,  2,  1,  2,  1,  2),
    source     = ""
  )
)

# Underlying databased on PERG SARS Systematic Review, published in the Lancet Microbe
# https://doi.org/10.1101/2024.08.13.24311934

datasets_SARS <- list(            # Severe Acute Respiratory Syndrome (SARS-CoV-1)

)


datasets_MERS <- list(            # Middle East Respiratory Syndrome (MERS-CoV)

)

# Underlying databased on PERG Zika Systematic Review, published in Nature Health
# https://doi.org/10.1038/s44360-025-00051-4

datasets_Zika <- list(            # Zika Virus Disease

)


datasets_Measles <- list(         # Measles
  d1  = list(
    # n defaults to sum(freq_count) = 116
    freq_value = c( 6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25),
    freq_count = c( 1,  5,  7,  9, 17, 18, 16,  6,  8, 10,  7,  1,  2,  4,  4,  1),
    source     = ""
  ),
  d2  = list(
    # n defaults to sum(freq_count) = 26
    freq_value = c( 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 25),
    freq_count = c( 1,  1,  2,  3,  3,  5,  1,  1,  3,  1,  1,  2,  1),
    source     = ""
  ),
  d3  = list(
    # n defaults to sum(freq_count) = 16
    freq_value = c(11, 12, 13, 14, 15, 16, 18, 22),
    freq_count = c( 1,  1,  3,  3,  4,  2,  1,  1),
    source     = ""
  ),
  d4  = list(median = 14, min = 12, max = 18, n =  8, source = ""),
  d5  = list(median = 12, min =  8, max = 17, n = 21, source = ""),
  d6  = list(median = 16, min = 11, max = 21, n = 30, source = ""),
  d7  = list(median = 14, min = 11, max = 22, n = 16, source = ""),
  d8  = list(median = 14, min = 10, max = 21, n = 34, source = ""),
  d9  = list(mean   = 13.8, sd = 2.7,         n = 22, source = ""),
  d10 = list(mean   = 14.2, sd = 2.9,         n = 38, source = ""),
  d11 = list(median = 14,   min = 13, max = 17, n =  9, source = "")   # QA concern: see PHAC dataset notes
)


datasets_Oropouche <- list(       # Oropouche Fever

)


datasets_Mpox <- list(            # Mpox

)


datasets_Cholera <- list(         # Cholera

)
