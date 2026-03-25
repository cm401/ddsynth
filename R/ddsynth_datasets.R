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
# Frequency tables have been extracted from studies with 10 or more patients only (as reported in the paper, though might be less as usuable incubation period data)

datasets_SARS <- list(            # Severe Acute Respiratory Syndrome (SARS-CoV-1)
  d1  = list(mean   =  5.9, sd  =  3.5,             n = 96,  source = "Wu (2003)"),
  d2  = list(mean   =  4.7, sd  =  4.6,             n = 234, source = "Virlogeux (2015), doi: doi.org/10.1097/ede.0000000000000339"),
  d3  = list(median =  4.0, min =  2.0, max = 10.0, n =  42, source = "Varia (2003)"),
  d4  = list(median =  3.0, min =  2.0, max =  6.0, n =  11, source = "Wong (2004), doi: doi.org/10.3201/eid1002.030452"),
  d5  = list(median =  5.0, min =  1.0, max = 15.0, n =   7, source = "Scales (2003), doi: doi.org/10.3201/eid0910.030525"),
  d6  = list(median =  4.0, min =  1.0, max = 18.0, n =  19, source = "Meltzer (2004), doi: doi.org/10.3201/eid1002.030426"),
  d7  = list(mean   =  5.3, sd  =  4.5,             n =  85, source = "McBryde (2006), doi: doi.org/10.1007/s11538-005-9005-4"),
  d8  = list(median =  6.0, min =  1.0, max = 15.0, n =  98, source = "Liu (2016), doi: doi.org/10.1371/journal.pone.0149988"),
  d9  = list(median =  6.0, min =  2.0, max = 16.0, n = 138, source = "Lee (2003), doi: doi.org/10.1056/NEJMoa030685"),
  d10 = list(median =  4.0, min =  2.0, max =  8.0, n =   7, source = "Hsu (2003), doi: doi.org/10.3201/eid0906.030264"),
  d11 = list(median =  7.0, min =  4.0, max = 12.0, n =  13, source = "Hsu (2003), doi: doi.org/10.3201/eid0906.030264"),
  d12 = list(mean   =  5.1, sd  =  2.2,             n =  50, source = "Goh (2006)"),
  d13 = list(median =  4.0, min =  3.0, max =  6.0, n =  32, source = "Chen (2003), doi: doi.org/10.1001/archotol.129.11.1157"),
  d14 = list(median =  6.0, Q1  =  3.0, Q3  = 10.0, n = 144, source = "Booth (2003), doi: doi.org/10.1001/jama.289.21.JOC30885"),
  d15 = list(mean   =  4.0, sd  =  3.0,             n =   4, source = "Avendano (2003)"),
  d16 = list(mean   =  3.5, sd  =  3.0,             n =  10, source = "Avendano (2003)"),
  d17 = list(
    freq_lower = c( 2,  1,  1,  1,  1,  3,  3,  1,  1,  2,  6,  2,  1,  5,  5,  1,  2,  13,  7),
    freq_upper = c(12,  4,  4, 11, 14,  3, 10,  6,  2,  2,  6,  6, 11, 11, 11,  5,  7,  18, 12),
    freq_count = c( 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1),
    source     = "Meltzer (2004), doi: 10.3201/eid1002.030426"
  ),
  d18 = list(
    # Note: four cases have freq_lower = 0 (unknown exposure start). Replaced with 0.1 to
    # avoid a NaN gradient in the Weibull likelihood (0^phi * log(0) = 0 * -Inf = NaN).
    # Lognormal and gamma are unaffected; this has negligible impact on all posteriors.
    freq_lower = c( 2.5,  4.5,  0.1,  0.5,  5.5,  2.5,  1.5,  3.5,  7.5,  1.5,  2.5,  7.5,  4.5,  4.5,  3.5,
                    0.5,  0.5,  2.5,  2.5,  0.1,  2.5,  0.1, 11.5,  0.5,  0.5,  1.5,  5.5,  3.5,  2.5, 12.5,
                    1.5,  8.5,  7.5,  1.5,  5.5,  4.5,  3.5,  1.5,  0.5,  0.5,  3.5,  7.5,  5.5,  2.5,  1.5,
                    1.5,  0.5,  1.5,  1.5,  1.5,  0.5,  8.5,  0.5,  0.1,  2.5,  5.5,  9.5,  2.5,  7.5,  4.5,
                    3.5,  1.5,  3.5,  3.5,  3.5,  4.5,  0.5),
    freq_upper = c( 6.5,  7.5,  2.5,  3.5,  6.5,  5.5,  4.5,  5.5, 11.5,  4.5,  5.5, 12.5,  5.5,  5.5,  6.5,
                    1.5,  4.5,  4.5,  4.5,  1.5,  4.5,  4.5, 13.5,  2.5,  1.5,  2.5, 10.5,  5.5,  4.5, 14.5,
                    6.5, 11.5,  9.5,  2.5,  7.5,  9.5,  7.5,  2.5,  4.5,  1.5,  4.5,  8.5,  6.5,  7.5,  2.5,
                    2.5,  1.5,  2.5,  2.5,  3.5,  1.5, 12.5,  1.5,  4.5,  4.5,  6.5, 10.5,  3.5,  8.5,  7.5,
                    5.5,  5.5,  8.5,  7.5,  7.5,  9.5,  3.5),
    freq_count = rep(1, 67),
    source     = "Farewell (2005), doi: 10.1002/sim.2206"
  ),
  d19 = list(
    freq_value = c( 3, 4, 5, 7, 8),
    freq_count = c( 5, 5, 3, 1, 1),
    source     = "Chow (2004), doi: 10.1136/bmj.37939.465729.44"
  ),
  d20 = list(
    freq_value = c( 2,  3, 4, 5, 6, 8),
    freq_count = c( 1, 10, 5 ,3, 2, 1),
    source     = "Olsen (2003), doi: 10.1056/NEJMoa031349"
  ),
  d21 = list(
    freq_value = c( 2, 3, 4, 5, 6),
    freq_count = c( 3, 3, 2, 2, 1),
    source     = "Wong (2004), doi: 10.3201/eid1002.030452"
  ),
  d22 = list(
    freq_lower = c( 4, 7, 6, 1, 3, 5, 10 ),
    freq_upper = c( 4, 7, 7, 5, 7, 6, 12 ),
    freq_count = c( 1, 1, 1, 1, 2, 1, 1),
    source     = "Dwosh (2003)"
  )
)

# Underlying databased on PERG MERS Systematic Review (currently unpublished)

datasets_MERS <- list(            # Middle East Respiratory Syndrome (MERS-CoV)
  d1  = list(median = 7.00, Q1  = 5.0, Q3  = 10.0, n =  73, source = "Cho (2016), doi: 10.1016/S0140-6736(16)30623-7"),
  d2  = list(median = 8.00, Q1  = 6.5, Q3  = 10.5, n =  14, source = "Nam (2017), doi: 10.1016/j.ijid.2017.02.008"),
  d3  = list(median = 4.00, Q1  = 3.0, Q3  =  8.0, n =  11, source = "Nam (2017), doi: 10.1016/j.ijid.2017.02.008"),
  d4  = list(mean   = 6.27, sd  = 4.35,            n =  18, source = "Al-Jasser (2019), doi: 10.1016/j.jiph.2018.09.008"), 
  d5  = list(median = 5.00, min = 2.0, max = 13.0, n =  36, source = "Kim (2015), doi: 10.1016/j.jiph.2018.09.008"),
  d6  = list(median = 7.00, min = 2.0, max = 14.0, n =  17, source = "Kim (2016), doi: 10.1016/j.phrp.2016.01.001"),
  d7  = list(median = 6.00, min = 2.0, max = 15.0, n =  36, source = "Park (2015), doi: 10.2807/1560-7917.es2015.20.25.21169"),
  d8  = list(median = 5.00, min = 2.0, max = 15.0, n =  92, source = "Liu (2016), doi: 10.1371/journal.pone.0149988"),
  d9 = list(median = 7.00, min = 3.0, max = 11.0, n = 128, source = "Liu (2016), doi: 10.1371/journal.pone.0149988"),
  d10 = list(      # Saudi Arabia
    freq_lower = c( 1, 1, 1, 1,  1, 2, 3, 4, 5, 6, 7, 8,  8, 9, 11, 14 ),
    freq_upper = c( 1, 4, 5, 6, 10, 2, 3, 8, 8, 6, 7, 8, 12, 9, 19, 14 ),
    freq_count = c( 1, 1, 2, 2,  1, 1, 6, 1, 1, 1, 1, 1,  1, 1,  1,  1 ),
    source     = "Assiri (2013), doi: 10.1056/NEJMoa1306742"
  ),
  d11 = list(      # South Korea
    freq_lower = c( 0.1, 0.1, 0.1, 0.1,  1,  2,  2,  2,  3,  3,  3,  4,  
                    4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  
                    6,  7,  7,  7,  7,  8,  8,  8,  8,  8,  9,  9, 10, 
                    10, 10, 11, 11, 12, 12, 12, 13, 15 ),
    freq_upper = c( 0.1,   6,  12,  15,  3,  3,  4,  8,  3,  5,  8,  4,  
                    5,  6,  8, 14,  5,  6,  7,  8,  9, 13,  6,  7,  8, 
                    12,  7,  8,  9, 11,  8,  9, 10, 11, 12,  9, 10, 10, 
                    12, 18, 12, 13, 14, 21, 27, 13, 17 ),
    freq_count = c(   1,   1,   1,   1,  4,  1,  3,  1,  5,  3,  1,  3,  
                      2,  4,  1,  1,  6,  4,  1,  2,  1,  1,  4,  1,  3,  
                      1,  2,  1,  3,  1,  6,  1,  2,  3,  1,  3,  2,  6,  
                      3,  1,  1,  1,  1,  1,  1,  1,  1 ),
    source     = "Virlogeux (2016), doi: 10.3201/eid2203.151437"
  )
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
  d11 = list(median = 14,   min = 13, max = 17, n =  9, source = "Sheline (1987)"),   # QA concern: see PHAC dataset notes
  d12 = list(
    freq_lower = c(12),
    freq_upper = c(15),
    freq_count = c(38),   # sum of the counts above
    source     = "Panum (1847), republished translation"
  )
)

# Underlying data based on mpox Systematic Review by Diaz Brochero et al (2025)
# published in BMJ Global Health: https://gh.bmj.com/content/10/1/e016906

datasets_Mpox <- list(            # Mpox
  d1  = list(median =  8, min =  2, max = 40, n =  78, source = "Angelo (2023), doi: 10.1016/S1473-3099(22)00651-X"),
  d2  = list(median =  6, Q1  =  4, Q3  =  9, n =  77, source = "Catala (2022), doi: 10.1111/bjd.21790"),
  d3  = list(median =  7, Q1  =  4, Q3  = 10, n = 179, source = "Choudury (2022), doi: 10.3238/arztebl.m2022.0340"),
  d4  = list(median =  9, min =  4, max = 15, n =  30, source = "Gaspari (2022), doi: 10.1128/jcm.01365-22"),
  d5  = list(mean   = 8.2, sd = 4.7,          n = 209, source = "Kroger (2023), doi: 10.1007/s15010-023-01997-x"),
  d6  = list(median =  6, Q1  =  3, Q3  =  8, n =  86, source = "Mailhe (2023), doi: 10.1016/j.cmi.2022.08.012"),
  d7  = list(mean   = 8.1, sd = 4.4,          n =  18, source = "Muira (2023), doi: 10.1093/infdis/jiad091"),
  d8  = list(median = 11, Q1  = 11, Q3  = 16, n =  16, source = "Moschese (2023), doi: 10.1016/j.jinf.2022.08.019"),
  d9  = list(median =  8, Q1  =  4, Q3  =  9, n =  18, source = "Nunez (2023), doi: 10.1016/j.lana.2022.100392"),
  d10 = list(median =  7, Q1  =  5, Q3  = 11, n = 181, source = "Tarin-Vicente (2022), doi: 10.1016/S0140-6736(22)01436-2"),
  d11 = list(median =  7, Q1  =  4, Q3  = 11, n =  51, source = "Thornhill (2022), doi: 10.1016/S0140-6736(22)02187-0"),
  d12 = list(median =  7, min =  3, max = 20, n =  23, source = "Thornhill (2022a), doi: 10.1056/NEJMoa2207323"),
  d13 = list(
    freq_lower = c( 4, 0.1,   1, 2,   3, 0.1,  9, 0.1, 10, 1,  6, 0.1,  7, 0.1,  8,  1,  2,  1, 0.1,  1,  3, 0.1 ),
    freq_upper = c( 6,   8,  23, 5,   6,  14, 11,  24, 12, 6,  8,  18, 16,   2, 17, 29, 10,  6,  23,  3,  8,  12 ),
    freq_count = c( 1,   1,   1, 1,   1,   1,  1,   1,  1, 1,  1,   1,  1,   1,  1,  1,  1,  1,   1,  1,  1,   1 ),
    source     = "Charniga (2022), doi: 10.1101/2022.06.22.22276713"
  ),
  d14 = list(
    freq_value = c(  2,  3,  5,  7,  8,  9, 10, 11, 17 ),
    freq_count = c(  2,  1,  4,  2,  2,  2,  4,  1,  1 ),
    source     = "Cobos (2023), doi: 10.37201/req/112.2022"
  ),
  d15 = list(
    freq_lower = c( 0.1, 0.1, 0.1, 0.1,  1,  2,  2,  3,  3,  4,  5,  6,  6,  6,  
                    6,  7,  7,  7,  8,  9,  9, 10, 11, 11, 12, 15, 19, 19 ),
    freq_upper = c(   2,   3,   7,  12,  4,  2,  9,  3,  8,  4,  5,  6, 12, 13, 
                      15,  7, 14, 14,  13,  9, 24,  10, 11, 11, 22, 22, 19, 24 ),
    freq_count = c(   1,   1,   1,   1,  1,  1,  1,  1,  1,  1,  1,  2,  1,  1,   
                      1,  2,  1,  1,   1,  2,  1,   1,  3,  1,  1,  1,  1,  1 ),
    source     = "Guzzetta (2022), doi: 10.3201/eid2810.221126"
  ),
  d16 = list(
    freq_lower = c( 0.1, 0.1, 0.1,   1,   1,   1,   2,   2,   2,   2,   2,   3,   3,   
                    3,   4,   4,   5,   5,   5,   6,   6,   6,   6,   7,   7,   7,   
                    7,   8,   8,   8,   9,   9,   9,   9,  10,  10,  10,  11,  11,  
                    12,  12,  13,  13,  14,  14,  15,  15,  15,  15,  16,  16,  17,  
                    18,  18,  18,  19,  19,  20,  20,  21,  21,  22,  30 ),
    freq_upper = c(   2,   3,   4,   1,   3,   6,   2,   3,   4,   6,   7,   3,   4,   
                      7,   4,   7,   5,   6,   8,   6,   7,  10,  11,   7,  10,  11,  
                      12,   8,   9,  13,   9,  10,  12,  14,  10,  11,  13,  11,  12,  
                      12,  15,  13,  15,  14,  19,  15,  16,  19,  20,  17,  20,  17,  
                      18,  19,  23,  19,  20,  21,  24,  23,  25,  26,  35 ),
    freq_count = c(   1,   1,   2,   1,   1,   1,   1,   1,   2,   1,   1,   2,   2,   
                      3,   7,   1,   6,   1,   1,   5,   2,   1,   1,   6,   1,   1,   
                      1,   4,   1,   4,   5,   1,   1,   1,   8,   2,   1,   4,   2,   
                      2,   3,   1,   1,   2,   2,   2,   1,   1,   4,   1,   1,   1,   
                      1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ),
    source     = "McFarland (2023), doi: 10.2807/1560-7917.ES.2023.28.27.2200806"
  )
)


datasets_Cholera <- list(         # Cholera

)
