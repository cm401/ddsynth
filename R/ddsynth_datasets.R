# =============================================================================
# ddsynth_datasets.R
# -----------------------------------------------------------------------------
# Built-in incubation period datasets for use with the ddsynth package.
#
# Each dataset contains summary statistics (central measure, variability, 
# and sample size) or frequency/interval-censored observations extracted from
# published epidemiological studies for a given pathogen. Data are organised
# as named lists (d1, d2, ...) within a per-pathogen list object.
#
# Supported summary types:
#   - Median + range (min/max)
#   - Median + IQR (Q1/Q3)
#   - Mean + standard deviation (sd)
#   - Frequency table (freq_value, freq_count)
#   - Interval-censored frequency table (freq_lower, freq_upper, freq_count)
#
# Underlying data are drawn from systematic reviews conducted by the
# Pathogen Epidemiology Review Group (PERG) at Imperial College London:
#   https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/related-initiatives/perg/ 
# for WHO priority pathogens, and from the broader literature for other pathogens, 
# utilising systematic reviews where available, and individual studies otherwise. 
#
# Each entry includes a 'source' field citing the original publication in the
# format: "Surname (year), doi: doi.org/DOI"
# =============================================================================

# Underlying data based on PERG Nipah Systematic Review (currently available as a pre-print), 
# https://doi.org/10.64898/2026.03.19.26348815

datasets_Nipah <- list(
  d1  = list(median = 10.0, min =  9.0, max = 12.0, n =  4, source = "Rahman (2012), doi: doi.org/10.1089/vbz.2011.0656"),
  d2  = list(median =  4.0, min =  2.0, max =  7.0, n =  6, source = "Rahman (2012), doi: doi.org/10.1089/vbz.2011.0656"),
  d3  = list(median =  9.0, min =  6.0, max = 14.0, n = 11, source = "Nikolay (2019), doi: doi.org/10.1056/NEJMoa1805376"),
  d4  = list(median =  9.5, min =  4.0, max = 14.0, n = 22, source = "Pallivalappil (2020), doi: doi.org/10.4103/jgid.jgid_4_19"),
  d5  = list(median =  8.0, min =  3.0, max = 20.0, n = 15, source = "Ching (2015), doi: doi.org/10.3201/eid2102.141433"),
  d6  = list(median =  9.0, min =  6.0, max = 11.0, n = 14, source = "Luby (2009), doi: doi.org/10.3201/eid1508.081237"),
  d7  = list(median = 10.0, min =  8.0, max = 15.0, n = 11, source = "Chandni (2019), doi: doi.org/10.1093/cid/ciz789"),
  d8  = list(median =  9.0, min =  6.0, max = 11.0, n = 11, source = "Hossain (2008), doi: doi.org/10.1086/529147"),
  d9  = list(median =  9.0, Q1  =  8.0, Q3  = 11.0, n = 82, source = "Nikolay (2019), doi: doi.org/10.1056/NEJMoa1805376"),
  d10 = list(mean   =  9.3, sd  =  1.9,             n = 18, source = "Thomas (2019), doi: doi.org/10.4103/ijcm.IJCM_198_19"),
  d11 = list(
    # n defaults to sum(freq_count) = 11
    freq_value = c(6, 7, 8, 9, 11, 12, 14),
    freq_count = c(1, 1, 3, 2,  1,  2,  1),
    source     = "Nikolay (2019), doi: doi.org/10.1056/NEJMoa1805376"
  )
)

# Underlying data based on PERG MVD Systematic Review, published in the Lancet ID
# https://doi.org/10.1016/S1473-3099(23)00515-7

datasets_MVD <- list(         # Marburg Virus Disease
  d1 = list(median =  7.0, min =  2.0, max = 14.0, n = 66, source = "Pavlin (2014), doi: doi.org/10.1186/1756-0500-7-906"),
  d2 = list(median = 10.0, Q1  =  8.0, Q3  = 13.0, n = 76, source = "Nsanzimana (2025), doi: doi.org/10.1056/NEJMoa2415816")
)

# Underlying data based on PERG EVD Systematic Review, published in the Lancet ID
# https://doi.org/10.1016/S1473-3099(24)00374-8

datasets_EVD <- list(         # Ebola Virus Disease
  d1  = list(median =  7.0, min =  2.0, max = 20.0, n =  116, source = "Wamala (2010), doi: doi.org/10.3201/eid1607.091525"),
  d2  = list(median =  6.0, min =  1.0, max = 16.0, n =   24, source = "Francesconi (2003), doi: doi.org/10.3201/eid0911.030339"),
  d3  = list(mean   =  9.2, sd  =  6.7,             n =   33, source = "Yan (2015), doi: doi.org/10.1007/s10096-015-2457-z"),
  d4  = list(mean   =  8.6, sd  =  6.1,             n =   20, source = "Yamin (2016), doi: doi.org/10.1016/j.ajic.2016.04.216"),
  d5  = list(mean   =  9.5, sd  =  4.0,             n =   76, source = "Muoghalu (2017), doi: doi.org/10.3389/fpubh.2017.00160"),
  d6  = list(mean   = 10.0, sd  =  1.0,             n =  291, source = "Lekone (2006), doi: doi.org/10.1111/j.1541-0420.2006.00609.x"),  # NOTE: unusually low SD — flag for pre_inference_checks()
  d7  = list(mean   =  9.9, sd  =  5.5,             n =  152, source = "Faye (2015), doi: doi.org/10.1016/S1473-3099(14)71075-8"),
  d8  = list(mean   =  9.7, sd  =  3.7,             n =    8, source = "Ajelli (2015), doi: doi.org/10.1186/s12916-015-0524-z"),
  d9  = list(mean   =  9.3, sd  =  1.9,             n =   20, source = "Chan (2020), doi: doi.org/10.1098/rsif.2020.0498"),
  d10 = list(        # West Africa 2014
    # n defaults to sum(freq_count) = 143
    freq_value = c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                   21, 22, 23, 25, 27, 30, 32, 35, 38, 42),
    freq_count = c( 1,  2,  4,  6,  8,  9, 10, 11, 12, 13,
                   11, 10,  8,  7,  5,  4,  3,  3,  2,  2,
                    2,  1,  1,  1,  1,  1,  1,  1,  1,  1),
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  ),
  d11 = list(        # Guinea
    # n defaults to sum(freq_count) = 58
    freq_value = c( 2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
                   12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 30, 38),
    freq_count = c( 1,  2,  3,  4,  4,  4,  4,  5,  5,  4,
                    4,  3,  3,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1),
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  ),
  d12 = list(        # Liberia
    # n defaults to sum(freq_count) = 52
    freq_value = c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 25, 27, 32, 35, 42),
    freq_count = c( 1,  1,  1,  2,  2,  3,  4,  4,  4,  5,
                    4,  4,  3,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1),
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  ),
  d13 = list(        # Sierra Leone
    # n defaults to sum(freq_count) = 30
    freq_value = c( 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
    freq_count = c( 1,  1,  2,  2,  2,  3,  3,  3,  3,  2,  2,  2,  1,  1,  1,  1),
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  )
)

# Underlying data based on PERG Lassa Fever Systematic Review, published in the Lancet Global Health
# https://doi.org/10.1016/S2214-109X(24)00379-6

datasets_Lassa <- list(           # Lassa Fever
  # Single dataset — tau not identifiable; predictions computed at mean(loc_d).
  # See prepare_stan_data_from_datasets() documentation.
  d1 = list(
    # n defaults to sum(freq_count) = 15
    freq_lower = c( 1,  2,  2,  3,  4,  6,  7,  8, 10, 11),
    freq_upper = c(12, 13, 17, 18, 19, 13, 22, 12, 10, 15),
    freq_count = c( 1,  1,  1,  2,  2,  2,  1,  2,  1,  2),
    source     = "Carey (1972), doi: 10.1016/0035-9203(72)90271-4"
  )
)

# Underlying data based on PERG SARS Systematic Review, published in the Lancet Microbe
# https://doi.org/10.1101/2024.08.13.24311934

datasets_SARS <- list(            # Severe Acute Respiratory Syndrome (SARS-CoV-1)

)

# Underlying databased on PERG MERS Systematic Review (currently unpublished)

datasets_MERS <- list(            # Middle East Respiratory Syndrome (MERS-CoV)

)

# Underlying data based on PERG Zika Systematic Review, published in Nature Health
# https://doi.org/10.1038/s44360-025-00051-4

datasets_Zika <- list(            # Zika Virus Disease
  d1  = list(median = 6, min = 2, max = 10, n =  111, source = "Sharma (2019), doi: 10.4103/INJMS.INJMS_65_19"),
  d2  = list(
    # n defaults to sum(freq_count) = 15
    freq_lower = c( 3,  3,  7,  1,  5,  4,  3,  3,  4,  1,  5,  0,  2,  1,  0,  0,  1,  0,  0,  5,  0,  0,  0,  0),
    freq_upper = c( 4, 12, 21, 18, 30, 29, 11, 18, 17, 185, 16, 10,  9, 24,  9, 30,  9, 77, 35, 28, 12, 16, 13, 20),
    freq_count = c( 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2),
    source     = "Lessler (2016), doi: 10.2471/BLT.16.174540"
  )
)

# Underlying data based on PHAC Measles Systematic Review, currently unpublished

datasets_Measles <- list(         # Measles
  d1  = list(
    # n defaults to sum(freq_count) = 116
    freq_value = c( 6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25),
    freq_count = c( 1,  5,  7,  9, 17, 18, 16,  6,  8, 10,  7,  1,  2,  4,  4,  1),
    source     = "Goodall (1931), pmc.ncbi.nlm.nih.gov/articles/PMC2313398/"   # historical BMJ article
  ),
  d2  = list(
    # n defaults to sum(freq_count) = 26
    freq_value = c( 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 25),
    freq_count = c( 1,  1,  2,  3,  3,  5,  1,  1,  3,  1,  1,  2,  1),
    source     = "Goodall (1931), pmc.ncbi.nlm.nih.gov/articles/PMC2313398/"   # historical BMJ article
  ),
  d3  = list(
    # n defaults to sum(freq_count) = 16
    freq_value = c(11, 12, 13, 14, 15, 16, 18, 22),
    freq_count = c( 1,  1,  3,  3,  4,  2,  1,  1),
    source     = "Ehresmann (1995), doi: doi.org/10.1093/infdis/171.3.679"
  ),
  d4  = list(median = 14, min = 12, max = 18, n =  8, source = "Bloch (1985)"),
  d5  = list(median = 12, min =  8, max = 17, n = 21, source = "Fielding (2005)"),
  d6  = list(median = 16, min = 11, max = 21, n = 30, source = "Kobayashi (2020), doi: doi.org/10.1016/j.vaccine.2020.05.067"),
  d7  = list(median = 14, min = 11, max = 22, n = 16, source = "Ehresmann (1995), doi: doi.org/10.1093/infdis/171.3.679"),
  d8  = list(median = 14, min = 10, max = 21, n = 34, source = "Papania (1999), doi: doi.org/10.1542/peds.104.5.e59"),
  d9  = list(mean   = 13.8, sd = 2.7,         n = 22, source = "Komabayashi (2018), doi: doi.org/10.7883/yoken.JJID.2018.083"),
  d10 = list(mean   = 14.2, sd = 2.9,         n = 38, source = "Komabayashi (2018), doi: doi.org/10.7883/yoken.JJID.2018.083"),
  d11 = list(median = 14,   min = 13, max = 17, n =  9, source = "Sheline (1987)")   # QA concern: see PHAC dataset notes
)


datasets_Oropouche <- list(       # Oropouche Fever

)


datasets_Mpox <- list(            # Mpox

)


datasets_Cholera <- list(         # Cholera

)
