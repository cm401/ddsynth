# ddsynth_datasets.R
# 
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
# 


# Nipah -------------------------------------------------------------------

# Underlying data based on PERG Nipah Systematic Review (currently available as a pre-print), 
# https://doi.org/10.64898/2026.03.19.26348815

datasets_Nipah <- list(
  d1  = list(median = 10.0, min =  9.0, max = 12.0, n =  4, subgroup = 'Bangladesh', country = 'Bangladesh',  source = "Rahman (2012), doi: 10.1089/vbz.2011.0656"),
  d2  = list(median =  4.0, min =  2.0, max =  7.0, n =  6, subgroup = 'Bangladesh', country = 'Bangladesh',  source = "Rahman (2012), doi: 10.1089/vbz.2011.0656"),
  d3  = list(median =  9.0, min =  6.0, max = 14.0, n = 11, subgroup = 'Bangladesh', country = 'Bangladesh',  source = "Nikolay (2019), doi: 10.1056/NEJMoa1805376"),
  d4  = list(median =  9.5, min =  4.0, max = 14.0, n = 22, country = 'India',       source = "Pallivalappil (2020), doi: 10.4103/jgid.jgid_4_19"),
  d5  = list(median =  8.0, min =  3.0, max = 20.0, n = 15, country = 'Philippines', source = "Ching (2015), doi: 10.3201/eid2102.141433"),
  d6  = list(median =  9.0, min =  6.0, max = 11.0, n = 14, subgroup = 'Bangladesh', country = 'Bangladesh',  source = "Luby (2009), doi: 10.3201/eid1508.081237"),
  d7  = list(median = 10.0, min =  8.0, max = 15.0, n = 11, country = 'India',       source = "Chandni (2019), doi: 10.1093/cid/ciz789"),
  d8  = list(median =  9.0, min =  6.0, max = 11.0, n = 11, subgroup = 'Bangladesh', country = 'Bangladesh',  source = "Hossain (2008), doi: 10.1086/529147"),
  d9  = list(median =  9.0, Q1  =  8.0, Q3  = 11.0, n = 82, subgroup = 'Bangladesh', country = 'Bangladesh',  source = "Nikolay (2019), doi: 10.1056/NEJMoa1805376"),
  d10 = list(mean   =  9.3, sd  =  1.9,             n = 18, country = 'India',       source = "Thomas (2019), doi: 10.4103/ijcm.IJCM_198_19"),
  d11 = list(
    # n defaults to sum(freq_count) = 11
    freq_value = c(6, 7, 8, 9, 11, 12, 14),
    freq_count = c(1, 1, 3, 2,  1,  2,  1),
    country = 'Bangladesh',
    subgroup = 'Bangladesh',
    source     = "Nikolay (2019), doi: 10.1056/NEJMoa1805376"
  )
)


# MVD ---------------------------------------------------------------------

# Underlying data based on PERG MVD Systematic Review, published in the Lancet ID
# https://doi.org/10.1016/S1473-3099(23)00515-7

datasets_MVD <- list(         # Marburg Virus Disease
  d1 = list(median =  7.0, min =  2.0, max = 14.0, n = 66, country = 'Mixed', source = "Pavlin (2014), doi: 10.1186/1756-0500-7-906"),
  d2 = list(median = 10.0, Q1  =  8.0, Q3  = 13.0, n = 76, country = 'Rwanda',   source = "Nsanzimana (2025), doi: 10.1056/NEJMoa2415816")
)


# EVD ---------------------------------------------------------------------

# Underlying data based on PERG EVD Systematic Review, published in the Lancet ID
# https://doi.org/10.1016/S1473-3099(24)00374-8

datasets_EVD <- list(         # Ebola Virus Disease
  d1  = list(median =  7.0, min =  2.0, max = 20.0, n =  116, country = 'Uganda',       subgroup = 'non-West Africa', source = "Wamala (2010), doi: 10.3201/eid1607.091525"),
  d2  = list(median =  6.0, min =  1.0, max = 16.0, n =   24, country = 'Uganda',       subgroup = 'non-West Africa', source = "Francesconi (2003), doi: 10.3201/eid0911.030339"),
  d3  = list(mean   =  9.2, sd  =  6.7,             n =   33, country = 'Sierra Leone', subgroup = 'West Africa',     source = "Yan (2015), doi: 10.1007/s10096-015-2457-z"),
  d4  = list(mean   =  8.6, sd  =  6.1,             n =   20, country = 'Sierra Leone', subgroup = 'West Africa',     source = "Yamin (2016), doi: 10.1016/j.ajic.2016.04.216"),
  d5  = list(mean   =  9.5, sd  =  4.0,             n =   76, country = 'Sierra Leone', subgroup = 'West Africa',     source = "Muoghalu (2017), doi: 10.3389/fpubh.2017.00160"),
  d6  = list(mean   = 10.0, sd  =  1.0,             n =  291, country = 'DRC',          subgroup = 'non-West Africa', source = "Lekone (2006), doi: 10.1111/j.1541-0420.2006.00609.x"),  # NOTE: unusually low SD — flag for pre_inference_checks()
  d7  = list(mean   =  9.9, sd  =  5.5,             n =  152, country = 'Guinea',       subgroup = 'West Africa',     source = "Faye (2015), doi: 10.1016/S1473-3099(14)71075-8"),
  d8  = list(mean   =  9.7, sd  =  3.7,             n =    8, country = 'Sierra Leone', subgroup = 'West Africa',     source = "Ajelli (2015), doi: 10.1186/s12916-015-0524-z"),
  d9  = list(mean   =  9.3, sd  =  1.9,             n =   20, country = 'Nigeria',      subgroup = 'non-West Africa', source = "Chan (2020), doi: 10.1098/rsif.2020.0498"),
  d10 = list(        # West Africa 2014
    # n defaults to sum(freq_count) = 143
    freq_value = c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                   21, 22, 23, 25, 27, 30, 32, 35, 38, 42),
    freq_count = c( 1,  2,  4,  6,  8,  9, 10, 11, 12, 13,
                   11, 10,  8,  7,  5,  4,  3,  3,  2,  2,
                    2,  1,  1,  1,  1,  1,  1,  1,  1,  1),
    country = 'Mixed',
    subgroup = 'West Africa', 
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  ),
  d11 = list(        # Guinea
    # n defaults to sum(freq_count) = 58
    freq_value = c( 2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
                   12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 30, 38),
    freq_count = c( 1,  2,  3,  4,  4,  4,  4,  5,  5,  4,
                    4,  3,  3,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1),
    country = 'Guinea',
    subgroup = 'West Africa', 
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  ),
  d12 = list(        # Liberia
    # n defaults to sum(freq_count) = 52
    freq_value = c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 25, 27, 32, 35, 42),
    freq_count = c( 1,  1,  1,  2,  2,  3,  4,  4,  4,  5,
                    4,  4,  3,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1),
    country = 'Liberia',
    subgroup = 'West Africa', 
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  ),
  d13 = list(        # Sierra Leone
    # n defaults to sum(freq_count) = 30
    freq_value = c( 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
    freq_count = c( 1,  1,  2,  2,  2,  3,  3,  3,  3,  2,  2,  2,  1,  1,  1,  1),
    country = 'Sierra Leone',
    subcountry = 'West Africa', 
    source     = "WHO Ebola Response Team (2014), doi: 10.1056/NEJMoa1411100"
  )
)


# Lassa -------------------------------------------------------------------

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
    country = "Nigeria",
    source     = "Carey (1972), doi: 10.1016/0035-9203(72)90271-4"
  )
)


# SARS --------------------------------------------------------------------

# Underlying data based on PERG SARS Systematic Review, published in the Lancet Microbe
# https://doi.org/10.1101/2024.08.13.24311934
# Frequency tables have been extracted from studies with 10 or more patients only (as reported in the paper, though might be less as usuable incubation period data)

datasets_SARS <- list(            # Severe Acute Respiratory Syndrome (SARS-CoV-1)
  d1  = list(mean   =  5.9, sd  =  3.5,             n = 96,  subgroup = "China", country = "China", source = "Wu (2003)"),
  d2  = list(mean   =  4.7, sd  =  4.6,             n = 234, subgroup = "Hong Kong", country = "Hong Kong", source = "Virlogeux (2015), doi: 10.1097/ede.0000000000000339"),
  d3  = list(median =  4.0, min =  2.0, max = 10.0, n =  42, subgroup = "Canada", country = "Canada", source = "Varia (2003)"),
  d4  = list(median =  3.0, min =  2.0, max =  6.0, n =  11, subgroup = "Hong Kong", country = "Hong Kong", source = "Wong (2004), doi: 10.3201/eid1002.030452"),
  d5  = list(median =  5.0, min =  1.0, max = 15.0, n =   7, subgroup = "Canada", country = "Canada", source = "Scales (2003), doi: 10.3201/eid0910.030525"),
  d6  = list(median =  4.0, min =  1.0, max = 18.0, n =  19, country = "Mixed", source = "Meltzer (2004), doi: 10.3201/eid1002.030426"),
  d7  = list(mean   =  5.3, sd  =  4.5,             n =  85, subgroup = "China",country = "China", source = "McBryde (2006), doi: 10.1007/s11538-005-9005-4"),
  d8  = list(median =  6.0, min =  1.0, max = 15.0, n =  98, subgroup = "Taiwan", country = "Taiwan", source = "Liu (2016), doi: 10.1371/journal.pone.0149988"),
  d9  = list(median =  6.0, min =  2.0, max = 16.0, n = 138, subgroup = "Hong Kong", country = "Hong Kong", source = "Lee (2003), doi: 10.1056/NEJMoa030685"),
  d10 = list(median =  4.0, min =  2.0, max =  8.0, n =   7, subgroup = "Singapore", country = "Singapore", source = "Hsu (2003), doi: 10.3201/eid0906.030264"),
  d11 = list(median =  7.0, min =  4.0, max = 12.0, n =  13, subgroup = "Singapore", country = "Singapore", source = "Hsu (2003), doi: 10.3201/eid0906.030264"),
  d12 = list(mean   =  5.1, sd  =  2.2,             n =  50, subgroup = "Singapore", country = "Singapore", source = "Goh (2006)"),
  d13 = list(median =  4.0, min =  3.0, max =  6.0, n =  32, subgroup = "Taiwan", country = "Taiwan", source = "Chen (2003), doi: 10.1001/archotol.129.11.1157"),
  d14 = list(median =  6.0, Q1  =  3.0, Q3  = 10.0, n = 144, subgroup = "Canada", country = "Canada", source = "Booth (2003), doi: 10.1001/jama.289.21.JOC30885"),
  d15 = list(mean   =  4.0, sd  =  3.0,             n =   4, subgroup = "Canada", country = "Canada", source = "Avendano (2003)"),
  d16 = list(mean   =  3.5, sd  =  3.0,             n =  10, subgroup = "Canada", country = "Canada", source = "Avendano (2003)"),
  d17 = list(
    freq_lower = c( 2,  1,  1,  1,  1,  3,  3,  1,  1,  2,  6,  2,  1,  5,  5,  1,  2,  13,  7),
    freq_upper = c(12,  4,  4, 11, 14,  3, 10,  6,  2,  2,  6,  6, 11, 11, 11,  5,  7,  18, 12),
    freq_count = c( 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,  1),
    country = "Mixed", 
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
    subgroup = "Hong Kong", 
    country = "Hong Kong", 
    source     = "Farewell (2005), doi: 10.1002/sim.2206"
  ),
  d19 = list(
    freq_value = c( 3, 4, 5, 7, 8),
    freq_count = c( 5, 5, 3, 1, 1),
    subgroup = "Singapore", 
    country = "Singapore", 
    source     = "Chow (2004), doi: 10.1136/bmj.37939.465729.44"
  ),
  d20 = list(
    freq_value = c( 2,  3, 4, 5, 6, 8),
    freq_count = c( 1, 10, 5 ,3, 2, 1),
    subgroup = "China",
    country = "China", 
    source     = "Olsen (2003), doi: 10.1056/NEJMoa031349"
  ),
  d21 = list(
    freq_value = c( 2, 3, 4, 5, 6),
    freq_count = c( 3, 3, 2, 2, 1),
    subgroup = "Hong Kong", 
    country = "Hong Kong", 
    source     = "Wong (2004), doi: 10.3201/eid1002.030452"
  ),
  d22 = list(
    freq_lower = c( 4, 7, 6, 1, 3, 5, 10 ),
    freq_upper = c( 4, 7, 7, 5, 7, 6, 12 ),
    freq_count = c( 1, 1, 1, 1, 2, 1, 1),
    subgroup = "Canada", 
    country = "Canada", 
    source     = "Dwosh (2003), PMID: 12771070"
  )
)


# MERS --------------------------------------------------------------------

# Underlying databased on PERG MERS Systematic Review (currently unpublished)

datasets_MERS <- list(            # Middle East Respiratory Syndrome (MERS-CoV)
  d1  = list(median = 7.00, Q1  = 5.0, Q3  = 10.0, n =  73, subgroup = 'Republic of Korea', country = 'Republic of Korea', source = "Cho (2016), doi: 10.1016/S0140-6736(16)30623-7"),
  d2  = list(median = 8.00, Q1  = 6.5, Q3  = 10.5, n =  14, subgroup = 'Republic of Korea', country = 'Republic of Korea', source = "Nam (2017), doi: 10.1016/j.ijid.2017.02.008"),
  d3  = list(median = 4.00, Q1  = 3.0, Q3  =  8.0, n =  11, subgroup = 'Republic of Korea', country = 'Republic of Korea', source = "Nam (2017), doi: 10.1016/j.ijid.2017.02.008"),
  d4  = list(mean   = 6.27, sd  = 4.35,            n =  18, subgroup = 'Saudi Arabia', country = 'Saudi Arabia',      source = "Al-Jasser (2019), doi: 10.1016/j.jiph.2018.09.008"), 
  d5  = list(median = 5.00, min = 2.0, max = 13.0, n =  36, subgroup = 'Republic of Korea', country = 'Republic of Korea', source = "Kim (2015), doi: 10.1177/1010539515610036"),
  d6  = list(median = 7.00, min = 2.0, max = 14.0, n =  17, subgroup = 'Republic of Korea', country = 'Republic of Korea', source = "Kim (2016), doi: 10.1016/j.phrp.2016.01.001"),
  d7  = list(median = 6.00, min = 2.0, max = 15.0, n =  36, subgroup = 'Republic of Korea', country = 'Republic of Korea', source = "Park (2015), doi: 10.2807/1560-7917.es2015.20.25.21169"),
  d8  = list(median = 5.00, min = 2.0, max = 15.0, n =  92, country = 'Middle East',       source = "Liu (2016), doi: 10.1371/journal.pone.0149988"),
  d9 = list(median = 7.00, min = 3.0, max = 11.0, n = 128,  country = 'Middle East',       source = "Liu (2016), doi: 10.1371/journal.pone.0149988"),
  d10 = list(      # Saudi Arabia
    freq_lower = c( 1, 1, 1, 1,  1, 2, 3, 4, 5, 6, 7, 8,  8, 9, 11, 14 ),
    freq_upper = c( 1, 4, 5, 6, 10, 2, 3, 8, 8, 6, 7, 8, 12, 9, 19, 14 ),
    freq_count = c( 1, 1, 2, 2,  1, 1, 6, 1, 1, 1, 1, 1,  1, 1,  1,  1 ),
    subgroup = 'Saudi Arabia', 
    country = 'Saudi Arabia', 
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
    subgroup = 'Republic of Korea', 
    country = 'Republic of Korea',
    source     = "Virlogeux (2016), doi: 10.3201/eid2203.151437"
  )
)


# Zika --------------------------------------------------------------------

# Underlying data based on PERG Zika Systematic Review, published in Nature Health
# https://doi.org/10.1038/s44360-025-00051-4

datasets_Zika <- list(            # Zika Virus Disease
  d1  = list(median = 6, min = 2, max = 10, n =  111, country = "India", source = "Sharma (2019), doi: 10.4103/INJMS.INJMS_65_19"),
  d2  = list(
    # n defaults to sum(freq_count) = 15
    freq_lower = c( 3,  3,  7,  1,  5,  4,  3,  3,  4,  1,  5,  0,  2,  1,  0,  0,  1,  0,  0,  5,  0,  0,  0,  0),
    freq_upper = c( 4, 12, 21, 18, 30, 29, 11, 18, 17, 185, 16, 10,  9, 24,  9, 30,  9, 77, 35, 28, 12, 16, 13, 20),
    freq_count = c( 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2),
    country    = "Mixed",
    source     = "Lessler (2016), doi: 10.2471/BLT.16.174540"
  )
)


# Measles -----------------------------------------------------------------

# Underlying data based on PHAC Measles Systematic Review, currently unpublished

datasets_Measles <- list(         # Measles
  d1  = list(
    # n defaults to sum(freq_count) = 116
    freq_value = c( 6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25),
    freq_count = c( 1,  5,  7,  9, 17, 18, 16,  6,  8, 10,  7,  1,  2,  4,  4,  1),
    country = 'United Kingdom',
    source     = "Goodall (1931), pmc.ncbi.nlm.nih.gov/articles/PMC2313398/"   # historical BMJ article
  ),
  d2  = list(
    # n defaults to sum(freq_count) = 26
    freq_value = c( 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 25),
    freq_count = c( 1,  1,  2,  3,  3,  5,  1,  1,  3,  1,  1,  2,  1),
    country = 'United Kingdom',
    source     = "Goodall (1931), pmc.ncbi.nlm.nih.gov/articles/PMC2313398/"   # historical BMJ article
  ),
  d3  = list(
    # n defaults to sum(freq_count) = 16
    freq_value = c(11, 12, 13, 14, 15, 16, 18, 22),
    freq_count = c( 1,  1,  3,  3,  4,  2,  1,  1),
    country = 'USA',
    source     = "Ehresmann (1995), doi: 10.1093/infdis/171.3.679"
  ),
  d4  = list(median = 14, min = 12, max = 18, n =  8, country = 'USA',       source = "Bloch (1985)"),
  d5  = list(median = 12, min =  8, max = 17, n = 21, country = 'Australia', source = "Fielding (2005)"),
  d6  = list(median = 16, min = 11, max = 21, n = 30, country = 'Japan',     source = "Kobayashi (2020), doi: 10.1016/j.vaccine.2020.05.067"),
  d7  = list(median = 14, min = 11, max = 22, n = 16, country = 'USA',       source = "Ehresmann (1995), doi: 10.1093/infdis/171.3.679"),
  d8  = list(median = 14, min = 10, max = 21, n = 34, country = 'USA',       source = "Papania (1999), doi: 10.1542/peds.104.5.e59"),
  d9  = list(mean   = 13.8, sd = 2.7,         n = 22, country = 'Japan',     source = "Komabayashi (2018), doi: 10.7883/yoken.JJID.2018.083"),
  d10 = list(mean   = 14.2, sd = 2.9,         n = 38, country = 'Japan',     source = "Komabayashi (2018), doi: 10.7883/yoken.JJID.2018.083"),
  d11 = list(median = 14,   min = 13, max = 17, n =  9, country = 'USA',     source = "Sheline (1987)"),   # QA concern: see PHAC dataset notes
  d12 = list(
    freq_lower = c(12),
    freq_upper = c(15),
    freq_count = c(38),   # sum of the counts above
    country = 'Faroe Islands',
    source     = "Panum (1847), republished translation"
  )
)


# Mpox --------------------------------------------------------------------

# Underlying data based on mpox Systematic Review by Diaz Brochero et al (2025)
# published in BMJ Global Health: https://gh.bmj.com/content/10/1/e016906

datasets_Mpox <- list(            # Mpox
  d1  = list(median =  8, min =  2, max = 40, n =  78, country = 'Mixed',       subgroup = 'IIb', source = "Angelo (2023), doi: 10.1016/S1473-3099(22)00651-X"),
  d2  = list(median =  6, Q1  =  4, Q3  =  9, n =  77, country = 'Spain',       subgroup = 'IIb',source = "Catala (2022), doi: 10.1111/bjd.21790"),
  d3  = list(median =  7, Q1  =  4, Q3  = 10, n = 179, country = 'Germany',     subgroup = 'IIb',source = "Choudury (2022), doi: 10.3238/arztebl.m2022.0340"),
  d4  = list(median =  9, min =  4, max = 15, n =  30, country = 'Italy',       subgroup = 'IIb',source = "Gaspari (2022), doi: 10.1128/jcm.01365-22"),
  d5  = list(mean   = 8.2, sd = 4.7,          n = 209, country = 'Germany',     subgroup = 'IIb',source = "Kroger (2023), doi: 10.1007/s15010-023-01997-x"),
  d6  = list(median =  6, Q1  =  3, Q3  =  8, n =  86, country = 'France',      subgroup = 'IIb',source = "Mailhe (2023), doi: 10.1016/j.cmi.2022.08.012"),
  d7  = list(mean   = 8.1, sd = 4.4,          n =  18, country = 'Netherlands', subgroup = 'IIb',source = "Muira (2023), doi: 10.1093/infdis/jiad091"),
  d8  = list(median = 11, Q1  = 11, Q3  = 16, n =  16, country = 'Italy',       subgroup = 'IIb',source = "Moschese (2023), doi: 10.1016/j.jinf.2022.08.019"),
  d9  = list(median =  8, Q1  =  4, Q3  =  9, n =  18, country = 'Mexico',      subgroup = 'IIb',source = "Nunez (2023), doi: 10.1016/j.lana.2022.100392"),
  d10 = list(median =  7, Q1  =  5, Q3  = 11, n = 181, country = 'Spain',       subgroup = 'IIb',source = "Tarin-Vicente (2022), doi: 10.1016/S0140-6736(22)01436-2"),
  d11 = list(median =  7, Q1  =  4, Q3  = 11, n =  51, country = 'Mixed',       subgroup = 'IIb',source = "Thornhill (2022), doi: 10.1016/S0140-6736(22)02187-0"),
  d12 = list(median =  7, min =  3, max = 20, n =  23, country = 'Mixed',       subgroup = 'IIb',source = "Thornhill (2022a), doi: 10.1056/NEJMoa2207323"),
  d13 = list(
    freq_lower = c( 4, 0.1,   1, 2,   3, 0.1,  9, 0.1, 10, 1,  6, 0.1,  7, 0.1,  8,  1,  2,  1, 0.1,  1,  3, 0.1 ),
    freq_upper = c( 6,   8,  23, 5,   6,  14, 11,  24, 12, 6,  8,  18, 16,   2, 17, 29, 10,  6,  23,  3,  8,  12 ),
    freq_count = c( 1,   1,   1, 1,   1,   1,  1,   1,  1, 1,  1,   1,  1,   1,  1,  1,  1,  1,   1,  1,  1,   1 ),
    subgroup = 'IIb',
    country = 'USA',
    source     = "Charniga (2022), doi: 10.1101/2022.06.22.22276713"
  ),
  d14 = list(
    freq_value = c(  2,  3,  5,  7,  8,  9, 10, 11, 17 ),
    freq_count = c(  2,  1,  4,  2,  2,  2,  4,  1,  1 ),
    subgroup = 'IIb',
    country = 'Spain',
    source     = "Cobos (2023), doi: 10.37201/req/112.2022"
  ),
  d15 = list(
    freq_lower = c( 0.1, 0.1, 0.1, 0.1,  1,  2,  2,  3,  3,  4,  5,  6,  6,  6,  
                    6,  7,  7,  7,  8,  9,  9, 10, 11, 11, 12, 15, 19, 19 ),
    freq_upper = c(   2,   3,   7,  12,  4,  2,  9,  3,  8,  4,  5,  6, 12, 13, 
                      15,  7, 14, 14,  13,  9, 24,  10, 11, 11, 22, 22, 19, 24 ),
    freq_count = c(   1,   1,   1,   1,  1,  1,  1,  1,  1,  1,  1,  2,  1,  1,   
                      1,  2,  1,  1,   1,  2,  1,   1,  3,  1,  1,  1,  1,  1 ),
    subgroup = 'IIb',
    country = 'Italy',
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
    subgroup = 'IIb',
    country = 'Mixed',
    source     = "McFarland (2023), doi: 10.2807/1560-7917.ES.2023.28.27.2200806"
  )
)


# Cholera -----------------------------------------------------------------

# Underlying data based on Cholera Systematic Review by Azman et al (2013)
# published in Journal of Infection: https://www.journalofinfection.com/article/S0163-4453(12)00347-7/fulltext

datasets_Cholera <- list(         # Cholera
  d1  = list(median =  21/24, min =  19/24, max = 49/24, n =  7,  country = 'USA', subgroup = "O1 Classical", source = "Morris (1995), doi: 10.1093/infdis/171.4.903"),
  d2  = list(median =  31/24, min =  19/24, max = 64/24, n =  11, country = 'USA', subgroup = "O1 Classical", source = "Morris (1995), doi: 10.1093/infdis/171.4.903"),
  d3  = list(
    freq_lower = c( 0.1, 0.1, 0.1, 0.1,   1, 1.0, 1.0, 1.5,   2,   2,   2,   3,   4,   8,  11 ),
    freq_upper = c( 1.0, 2.0, 3.0, 1.0, 1.0, 1.5, 3.0, 2.0,   2,   3,   2,   3,   4,   8,  11 ),
    freq_count = c(    6,   2,   2,   1,   1,   1,   1,   2,   2,   2,   1,   1,   1,   1,   1 ),
    subgroup = "O1 Classical",
    country = 'United Kingdom',
    source = "Snow (1854), 'On the Mode of Communication of Cholera'"
  ),
  d4  = list(
    # Histogram bins read in hours, converted to days (/24)
    # Bin at x = onset between x and x+10 hours
    freq_lower = c( 20,  30,  40,  50,  60,  70,  80,  90, 100, 110, 120, 150 ) / 24,  # days: 0.833, 1.25, 1.667, 2.083, 2.5, 2.917, 3.333, 3.75, 4.167, 4.583, 5.0, 6.25
    freq_upper = c( 30,  40,  50,  60,  70,  80,  90, 100, 110, 120, 130, 160 ) / 24,  # days: 1.25, 1.667, 2.083, 2.5, 2.917, 3.333, 3.75, 4.167, 4.583, 5.0, 5.417, 6.667
    freq_count = c(  3,   5,   5,   7,   4,   4,   3,   3,   2,   2,   1,   1 ),
    subgroup = "O1 Classical",
    country = 'USA',
    source = "Cash (1974), doi: 10.1093/infdis/129.1.45"
  ),
  d5  = list(
    # Histogram bins read in hours, converted to days (/24)
    # Bin at x = onset between x and x+10 hours
    freq_lower = c( 10,  20,  30,  40,  50,  60,  70,  80 ) / 24,  # days: 0.417, 0.833, 1.25, 1.667, 2.083, 2.5, 2.917, 3.333
    freq_upper = c( 20,  30,  40,  50,  60,  70,  80,  90 ) / 24,  # days: 0.833, 1.25, 1.667, 2.083, 2.5, 2.917, 3.333, 3.75
    freq_count = c(  3,   5,   1,   2,   7,   2,   1,   1 ),
    subgroup = "O1 Classical",
    country = 'USA',
    source = "Cash (1974), doi: 10.1093/infdis/129.1.45"
  ),
  d6  = list(
    freq_lower = c( 0.1,  1,  2,  3,  4,  5,  6 ),
    freq_upper = c(   1,  2,  3,  4,  5,  6,  7 ),
    freq_count = c(   6, 26, 23,  6,  3,  1,  1 ),
    country = 'Mixed',
    source = "Eberhart-Phillips (1996), doi: 10.1017/S0950268800058891"
  ),
  d7  = list(
    # Exact incubation period values in hours, converted to days (/24)
    freq_value = c(  4,  24,  31,  34,  36,  38,  44,  57,  59,  61,  62, 114, 133, 203 ) / 24,
    freq_count = c(  1,   1,   1,   3,   2,   2,   1,   1,   1,   1,   1,   1,   1,   1 ),
    subgroup = "O1 El Tor",
    country = 'Singapore',
    source = "Goh (1984), doi: 10.1093/ije/13.2.210"
  ),
  d8  = list(
    # Vaccinees who developed clinical cholera, ranges in hours converted to days (/24)
    # Table IV: Inaba 569B, 10^6 CFU — vaccinees (3/6 attacked): range 75-97 hrs
    # Table V:  Ogawa 395, 10^6 CFU  — vaccinees (6/8 attacked): range 15-110 hrs
    freq_lower = c( 75,  15 ) / 24,
    freq_upper = c( 97, 110 ) / 24,
    freq_count = c(  3,   6 ),
    subgroup = "O1 Classical",
    country = 'USA',
    source = "Levine (1979), doi: 10.1016/0035-9203(79)90119-6"
  ),
  d9  = list(
    # Unimmunized controls who developed clinical cholera, ranges in hours converted to days (/24)
    # Table IV: Inaba 569B, 10^6 CFU — controls (3/6 attacked): range 25-70 hrs
    # Table V:  Ogawa 395, 10^6 CFU  — controls (8/8 attacked): range 15-104 hrs
    freq_lower = c( 25,  15 ) / 24,
    freq_upper = c( 70, 104 ) / 24,
    freq_count = c(  3,   8 ),
    country = 'USA',
    subgroup = "O1 Classical",
    source = "Levine (1979), doi: 10.1016/0035-9203(79)90119-6"  
  ),
  d10 = list(
    # Exact dot plot values in hours, converted to days (/24)
    # All four gastric condition groups combined (circles, triangles, squares)
    # Gastrectomized: 0(x4), 2(x2), 5, 12(x2), 15, 22(x2), 24(x2), 48(x6), 72
    # Achlorhydric:   12(x2), 22(x4), 48(x2), 72, 96
    # Normochlorhydric: 22, 48, 72(x3), 96(x2), 144(x3)
    # Not examined:   22(x2), 48(x4), 72(x2), 96, 120
    freq_value = c( 0.1,  2,  5, 12, 15, 22, 24, 48, 72, 96, 120, 144 ) / 24,
    freq_count = c(   4,  2,  1,  4,  1, 11,  2, 13,   7,  4,   1,   3 ),
    country = 'Italy',
    subgroup = "O1 El Tor",
    source = "Schiraldi (1974), PMCID: PMC2366304"
  ),
  d11 = list(
    freq_lower = c(0.1, 24, 49, 72)/24,
    freq_upper = c(24, 48, 72, 150)/24,
    freq_count = c(1, 16, 6, 2),
    subgroup = "O1 El Tor",
    country = 'USA',
    source = "Sutton (1974), doi: 10.1017/S0022172400023688"
  )
)


# RVF ---------------------------------------------------------------------

# Underlying data based on VBD Systematic Review by Rudolph et al (2014)
# published in AJTMH: https://doi.org/10.4269/ajtmh.13-0403

datasets_RVF <- list(
    # 6 cases: males aged 32-64, South Africa (Johannesburg), contact with infected animal tissue during 1950-51 SA epizootic
    # Mundel B, Gear J. S Afr Med J. 1951;25:797-800. PMID:1488672
    # 1 case: adult male, USA, laboratory infection while grinding mouse brain tissue
    # Became ill 6 days after accidental aerosol/skin contamination
    # Sabin AB, Blumberg RW. Proc Soc Exp Biol Med. 1947;64:385-389. doi:10.3181/00379727-64-15803
    # Note: incubation of 6 days confirmed in secondary literature (PMC10316118)
    # 6 cases: NAMRU-3 staff, Egypt, November 1977
    # All 6 developed virologically confirmed RVF 3 days after
    # inhaling virus aerosol during sheep slaughter (single shared exposure)
    # Hoogstraal H et al. Trans R Soc Trop Med Hyg. 1979;73(6):624-629.
    # doi:10.1016/0035-9203(79)90005-1  PMID:44038
  d1 = list(
    freq_value = c(3,4,6),
    freq_count = c(6,5,1),
    country = 'Mixed',
    source = "Hoogstraal (1979), Sabin (1947), and Mundel & Gear (1951)"
  ),
  d2 = list(
    # Veterinary school cluster, South Africa 2008 RVF outbreak
    # Exposure date known (animal autopsy); mean = 4.3 days, range 4-5 days
    # Median estimated via inverse Wan et al. (2014) Scenario S1 formula:
    # m = [4n*xbar - (n+1)(a+b)] / [2(n-1)]; result ~4.0 days across n=5-7
    # min and max taken directly from paper
    median = 4.0, min    = 4, max    = 5, n      = 6,   # approximate; veterinary school subgroup of 8 confirmed cases
    country = 'South Africa',
    source = "Archer (2011), doi:10.7196/samj.4544"
  )
)

# CCHF --------------------------------------------------------------------

# PubMed search: (CCHF OR "Crimean-Congo hemorrhagic fever" ) AND "incubation period" 

datasets_CCHF <- list( 
  d1 <- list( mean   = 4.00, sd     = 2.40, n      = 49, source = "Arslan (2024), doi: 10.14744/nci.2023.09815" ), # Tick-bite patients only (55.7% of 88 total). Range 1-11d.
  d2 <- list(           #N=12 Iranian HCW, blood/nosocomial exposure only.
    freq_value  = c(22, 2, 2, 7, 7, 5, 13, 6, 4, 1, 5, 8),
    freq_count  = c( 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1),
    country     = 'Iran',
    source      = "Fazlalipour (2024), doi:10.1186/s12879-024-10199"
  ),
  d3 <- list(           #Table 2 footnotes. Nosocomial cases: Case A intact skin/blood IP<20h (encoded as 1d, anomalous per authors), Case B conjunctival splash IP=5d, Case C needlestick IP=8d.
    freq_value = c(1, 5, 8),
    freq_count = c(1, 1, 1),
    country    = 'Iran',
    source = "Naderi (2013), doi:10.4269/ajtmh.2012.12-0337"
  ),
  d4 <- list( median = 3, Q1 = 2, Q3 = 4, n = 64, country = 'Turkey', source = "Beştepe Dursun (2021), doi: 10.1002/jca.21875"),
  d5 <- list( mean   = 5.5, sd     = 3.60, n      = 64, country = 'Turkey', source = "Koksal (2010), doi: 10.1016/j.jcv.2009.11.007" ), 
  d6 <- list( mean   = 4.9, sd     = 3.90, n      = 72, country = 'Turkey', source = "Koksal (2010), doi: 10.1016/j.jcv.2009.11.007" )
)

datasets_CCHF_extended <- list(  # this includes approximation and one dataset purely of outliers.
  d1 <- list( mean   = 4.00, sd     = 2.40, n      = 49, 
              source = "Arslan (2024), doi:10.14744/nci.2023.09815",
              country = 'Turkey',
              subgroup = "tick-bite" ), # Tick-bite patients only (55.7% of 88 total). Range 1-11d.
  d2 <- list(           #N=12 Iranian HCW, blood/nosocomial exposure only.
    freq_value  = c(22, 2, 2, 7, 7, 5, 13, 6, 4, 1, 5, 8),
    freq_count  = c( 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1),
    source      = "Fazlalipour (2024), doi:10.1186/s12879-024-10199",
    country     = 'Iran',
    subgroup    = "nosocomial"
  ),
  d3 <- list(           #N=12 outlier tick-bite cases (IP>12d) selected from 312 tick-bite patients. OUTLIER SERIES 
    freq_value  = c(13, 14, 15, 15, 16, 20, 22, 23, 24, 27, 41, 53),
    freq_count  = c( 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1),
    country = 'Turkey',
    source      = "Kaya (2011), doi:10.1016/j.ijid.2011.03.007",
    subgroup    = "tick-bite"
  ),
  d4 <- list(           #Table 2 footnotes. Nosocomial cases: Case A intact skin/blood IP<20h (encoded as 1d, anomalous per authors), Case B conjunctival splash IP=5d, Case C needlestick IP=8d.
    freq_value = c(1, 5, 8),
    freq_count = c(1, 1, 1),
    country    = 'Iran',
    source = "Naderi (2013), doi:10.4269/ajtmh.2012.12-0337",
    subgroup = "nosocomial"
  ),
  d5 <- list( median = 3, Q1 = 2, Q3 = 4, n = 64, 
              source = "Beştepe Dursun (2021), doi: 10.1002/jca.21875",
              country = 'Turkey',
              subgroup = "other"),  # mix of patients in hospital
  d6 <- list( mean   = 5.5, sd     = 3.60, n      = 64, 
              country = 'Turkey',
              source = "Koksal (2010), doi: 10.1016/j.jcv.2009.11.007",
              subgroup = "other"),  # mix of patients in hospital, ribavirin + supportive care
  d7 <- list( mean   = 4.9, sd     = 3.90, n      = 72, 
              country = 'Turkey',
              source = "Koksal (2010), doi: 10.1016/j.jcv.2009.11.007",
              subgroup = "other"),  # mix of patients in hospital, supportive care only
  # --- Swanepoel et al. 1989 (Rev Infect Dis 11 Suppl 4:S794-S800, PMID 2749111)
  # 50 cases, South Africa 1981-1987 (superset of 1987 paper)
  cchf_swanepoel_1989_tick <- list(
    median = 2.5,   # Wan est. = 2.3, rounded up; mean 3.4, right-skewed [Median estimated from mean (3.4d) via Wan et al. (2014) BMC Med Res Methodol 14:135.]
    min    = 2.0,
    max    = 7.0,
    n      = 17,
    country = 'South Africa',
    source = "Swanepoel (1989), doi:10.1093/clinids/11.Supplement_4.S794",
    subgroup = "tick-bite"
  ),
  cchf_swanepoel_1989_livestock <- list(
    median = 4.0,   # Wan est. = 3.9
    min    = 2.0,
    max    = 9.0,
    n      = 14,
    country = 'South Africa',
    source = "Swanepoel (1989), doi:10.1093/clinids/11.Supplement_4.S794",
    subgroup = "other"  # livestock
  ),
  cchf_swanepoel_1989_nosocomial <- list(
    median = 6.0,   # Wan est. = 6.2, rounded to 6.0
    min    = 3.0,
    max    = 7.0,
    n      = 5,
    country = 'South Africa',
    source = "Swanepoel (1989), doi:10.1093/clinids/11.Supplement_4.S794",
    subgroup = "nosocomial"
  ),
  cchf_pshenichnaya_2016 <- list( # Table 1 and case narratives. N=3 secondary cases via possible sexual contact, Russia. IPs interval-censored from exposure window. Transmission route unconfirmed — extended dataset only.
    freq_lower = c(3, 2, 8),
    freq_upper = c(5, 3, 9),
    freq_count = c(1, 1, 1),
    country = 'Russia',
    source = "Pshenichnaya (2016), doi:10.1016/j.ijid.2016.02.1008",
    subgroup = "sexual-transmission"
  )
)


# COVID-19 ----------------------------------------------------------------

# Underlying data based on COVID-19 Systematic Review by Xu et al (2023)
# published in BMC Medicine: https://link.springer.com/article/10.1186/s12916-023-03070-8

datasets_COVID_19 <- list(
  
  # --- Wildtype ---
  d1  = list( median = 7,    min = 4,    max = 12,    n = 8,    source = "Shen (2020), doi: 10.1093/ofid/ofaa231",                        country = "China",       subgroup = "Wildtype" ),
  d2  = list( median = 4.3,  Q1 = 2.5,  Q3 = 6.5,    n = 53,   source = "Bender (2021), doi: 10.3201/eid2704.204576",                     country = "Germany",     subgroup = "Wildtype" ),
  d3  = list( median = 6,    min = 1,    max = 13,    n = 27,   source = "Liu (2020), doi: 10.1097/JCMA.0000000000000411",                country = "China",       subgroup = "Wildtype" ),
  d4  = list( mean = 8.23,   sd = 3.58,              n = 22,   source = "Song (2020), doi: 10.1016/j.jinf.2020.04.018",                   country = "China",       subgroup = "Wildtype" ),
  d5  = list( median = 3,    min = 0.1,  max = 15,    n = 10,   source = "Ki (2020), doi: 10.4178/epih.e2020007",                         country = "South Korea", subgroup = "Wildtype" ),
  d6  = list( median = 8.5,  min = 1,    max = 24,    n = 28,   source = "Mao (2020), doi: 10.1186/s12889-020-09606-4",                   country = "China",       subgroup = "Wildtype" ),
  d7  = list( median = 7,    min = 2,    max = 12,    n = 8,    source = "Zhang (2020), doi: 10.1186/s12879-020-05570-x",                 country = "China",       subgroup = "Wildtype" ),
  d8  = list( median = 8,    min = 4,    max = 13,    n = 23,   source = "Zhang (2020), doi: 10.1186/s12879-020-05570-x",                 country = "China",       subgroup = "Wildtype" ),
  d9  = list( median = 10,   min = 7,    max = 15,    n = 46,   source = "Zhang (2020), doi: 10.1186/s12879-020-05570-x",                 country = "China",       subgroup = "Wildtype" ),
  d10 = list( median = 9,    min = 6,    max = 13,    n = 85,   source = "Guo (2020), doi: 10.1186/s12916-020-01719-2",                   country = "China",       subgroup = "Wildtype" ),
  d11 = list( median = 5,    min = 2,    max = 8,     n = 2907, source = "Nie (2020), doi: 10.1093/infdis/jiaa211",                       country = "China",       subgroup = "Wildtype" ),
  d12 = list( median = 8.5,  min = 6,    max = 12,    n = 75,   source = "Du (2021), doi: 10.7883/yoken.JJID.2021.274",                   country = "China",       subgroup = "Wildtype" ),
  d13 = list( mean = 9.1,    sd = 3.7,               n = 43,   source = "Hua (2020), doi: 10.1002/jmv.26180",                             country = "China",       subgroup = "Wildtype" ),
  d14 = list( median = 5,    min = 1,    max = 11,    n = 15,   source = "Wong (2020), doi: 10.4269/ajtmh.20-0771",                       country = "Brunei",      subgroup = "Wildtype" ),
  d15 = list( median = 5.7,  Q1 = 3.2,  Q3 = 8.8,    n = 268,  source = "Hu (2021), doi: 10.1038/s41467-021-21710-6",                     country = "China",       subgroup = "Wildtype" ),
  d16 = list( mean = 6.8,    sd = 4.1,               n = 254,  source = "Zhao (2021a), doi: 10.1016/j.epidem.2021.100482",                country = "China",       subgroup = "Wildtype" ),
  
  d17 = list(
    # n=256, Bavaria Germany, Wildtype. Jan 20 – Mar 19 2020.
    # Interval-censored: lower = onset - right_expo, upper = onset - left_expo.
    # 67 cases where left_expo = right_expo (point exposure): lower set to 0.1.
    # 10 cases with upper > 20 days reflect wide exposure windows (uninformative but valid).
    # Replaces summary-statistic entry d17 (mean=4.6, sd=3.0, n=256).
    freq_lower = c(
      6, 3, 4, 6, 1, 2, 2, 2, 0.1, 3, 9, 0.1, 2, 0.1, 0.1, 2, 3, 4, 1, 1,
      1, 2, 1, 5, 0.1, 1, 7, 0.1, 0.1, 7, 0.1, 6, 1, 3, 2, 5, 0.1, 1, 2, 2,
      4, 0.1, 0.1, 12, 1, 0.1, 11, 9, 5, 3, 3, 7, 7, 7, 6, 6, 0.1, 1, 1, 1,
      3, 1, 1, 9, 2, 0.1, 0.1, 4, 2, 2, 4, 4, 1, 5, 3, 8, 15, 2, 2, 3,
      15, 3, 1, 2, 4, 4, 7, 0.1, 0.1, 3, 6, 1, 0.1, 0.1, 2, 1, 1, 1, 3, 3,
      0.1, 1, 3, 9, 2, 0.1, 2, 0.1, 7, 9, 5, 4, 0.1, 0.1, 0.1, 0.1, 3, 3, 1, 1,
      11, 2, 0.1, 2, 2, 1, 5, 6, 9, 0.1, 5, 0.1, 2, 1, 2, 7, 1, 1, 7, 0.1,
      5, 3, 1, 4, 0.1, 1, 3, 2, 1, 0.1, 0.1, 4, 1, 10, 0.1, 2, 2, 0.1, 5, 10,
      0.1, 1, 2, 2, 0.1, 0.1, 7, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 3, 6, 0.1, 0.1, 0.1, 0.1,
      2, 1, 2, 2, 1, 3, 0.1, 2, 1, 1, 1, 0.1, 2, 2, 2, 0.1, 1, 0.1, 3, 2,
      2, 0.1, 2, 2, 8, 1, 0.1, 1, 1, 4, 5, 1, 8, 15, 2, 0.1, 3, 2, 0.1, 4,
      2, 0.1, 0.1, 2, 0.1, 5, 1, 6, 0.1, 2, 1, 0.1, 1, 0.1, 0.1, 0.1, 2, 1, 2, 2,
      2, 2, 0.1, 1, 2, 7, 0.1, 0.1, 3, 0.1, 9, 9, 1, 6, 5, 2
    ),
    freq_upper = c(
      13, 6, 6, 10, 3, 4, 6, 43, 7, 10, 12, 3, 6, 2, 3, 5, 5, 7, 3, 8,
      4, 5, 4, 12, 5, 2, 13, 6, 7, 9, 2, 13, 8, 10, 7, 6, 3, 6, 8, 9,
      6, 26, 7, 12, 3, 11, 16, 15, 25, 6, 9, 11, 11, 9, 7, 10, 7, 5, 7, 3,
      9, 8, 3, 16, 6, 3, 8, 11, 5, 6, 4, 6, 8, 9, 10, 12, 17, 8, 9, 9,
      20, 6, 4, 8, 7, 8, 13, 6, 1, 9, 13, 3, 6, 7, 3, 3, 2, 23, 6, 10,
      5, 3, 7, 16, 6, 5, 3, 4, 8, 14, 10, 12, 7, 21, 3, 5, 9, 7, 6, 5,
      15, 8, 4, 9, 5, 8, 7, 9, 14, 4, 12, 4, 5, 5, 5, 8, 9, 6, 14, 7,
      5, 7, 3, 5, 5, 7, 8, 9, 5, 8, 6, 7, 9, 13, 2, 5, 6, 3, 6, 12,
      3, 4, 4, 9, 2, 5, 14, 4, 4, 5, 3, 5, 3, 6, 7, 8, 6, 4, 4, 5,
      5, 12, 5, 9, 5, 8, 6, 3, 6, 5, 3, 15, 4, 6, 22, 12, 4, 4, 32, 5,
      6, 36, 9, 5, 12, 4, 41, 3, 3, 7, 13, 4, 13, 17, 6, 4, 8, 5, 3, 11,
      6, 8, 7, 5, 2, 9, 5, 10, 3, 6, 3, 6, 7, 6, 4, 41, 4, 4, 4, 6,
      10, 7, 4, 5, 6, 7, 4, 5, 7, 19, 10, 13, 3, 8, 7, 4
    ),
    freq_count = c(rep(1, 256)),
    source  = "Böhm (2021), doi: 10.1017/S0950268821000510",
    country = "Germany",
    subgroup = "Wildtype"
  ),
  
  # d18 removed: superseded by d_bohmer_2020 below (same study, interval-censored individual data)
  d19 = list( median = 5.4,  min = 1,    max = 21,    n = 178,  source = "Yang (2020), doi: 10.1017/S0950268820001338",                   country = "China",       subgroup = "Wildtype" ),
  d20 = list( median = 5.6,  min = 1.35, max = 13.04, n = 19,   source = "Bui (2020), doi: 10.1371/journal.pone.0243889",                 country = "Vietnam",     subgroup = "Wildtype" ),
  d21 = list( median = 5.5,  min = 3,    max = 12,    n = 8,    source = "Han (2020), doi: 10.4178/epih.e2020056",                        country = "South Korea", subgroup = "Wildtype" ),
  
  d22 = list(
    # Fig 1, Kong (2020). Exact integer incubation periods, n=136.
    # All three age groups (0-14, 15-64, ≥65) combined.
    # Point observations: lower = upper = value; 0.1 floor not needed (min = 1).
    # Replaces summary-statistic entry d22 (median=8.3, Q1=4.9, Q3=12, n=136).
    freq_value = c(1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17),
    freq_count = c(2,  8,  6, 11,  8, 12, 13,  8, 11, 13,  7, 10,  9,  3,  5,  5,  5),
    source  = "Kong (2020), doi: 10.1002/agm2.12114",
    country = "China",
    subgroup = "Wildtype"
  ),
  
  d23 = list(
    # Fig 1 histogram, Xiao (2021). Shenzhen (n=176) + Hefei (n=41) combined, n=217.
    # Digitised from histogram; n=217 and mean=8.59 verified against reported 8.58.
    freq_value = c( 1,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18),
    freq_count = c( 4, 15, 19, 16, 26, 20, 15, 12, 21, 17, 13,  9,  8,  6, 12,  2,  2),
    source  = "Xiao (2021), doi: 10.1186/s12199-021-00935-3",
    country = "China",
    subgroup = "Wildtype"
  ),
  
  d24 = list( mean = 5.83,   sd = 3.73,              n = 180,  source = "Dai (2020), doi: 10.2147/RMHP.S257907",                          country = "China",       subgroup = "Wildtype" ),
  d25 = list( median = 2,    min = 1,    max = 4,     n = 6,    source = "Huang (2020), doi: 10.1016/j.jinf.2020.03.006",                 country = "China",       subgroup = "Wildtype" ),
  d26 = list( mean = 6.5,    sd = 4.58,              n = 102,  source = "Zhao (2021b), doi: 10.1097/MD.0000000000027846",                 country = "China",       subgroup = "Wildtype" ),
  d27 = list( mean = 8.67,   sd = 5.16,              n = 957,  source = "Li (2020), doi: 10.3389/fpubh.2020.577431",                      country = "China",       subgroup = "Wildtype" ),
  d28 = list( median = 4,    min = 0.1,  max = 11,    n = 45,   source = "Bernal Lopez (2022), doi: 10.2807/1560-7917.ES.2022.27.15.2001551",   country = "UK",          subgroup = "Wildtype" ),  # min raised 0→0.5: integer-day same-day onset recorded as 0
  d29 = list( median = 4,    min = 2,    max = 11,    n = 12,   source = "Bernal Lopez (2022), doi: 10.2807/1560-7917.ES.2022.27.15.2001551",   country = "UK",          subgroup = "Wildtype" ),
  
  # --- Omicron ---
  d30 = list( mean = 4.58,   sd = 1.72,              n = 57,   source = "Mefsin (2022), doi: 10.3201/eid2809.220613",                     country = "China",       subgroup = "Omicron" ),
  d31 = list( mean = 4.42,   sd = 1.42,              n = 23,   source = "Mefsin (2022), doi: 10.3201/eid2809.220613",                     country = "China",       subgroup = "Omicron" ),
  d32 = list( median = 3,    min = 0.1,  max = 8,     n = 81,   source = "Brandal (2021), doi: 10.2807/1560-7917.ES.2021.26.50.2101147",  country = "Norway",      subgroup = "Omicron" ),  # min raised 0→0.5: integer-day same-day onset recorded as 0
  d34 = list( mean = 3.2,    sd = 2.2,               n = 258,  source = "Backer (2022), doi: 10.2807/1560-7917.ES.2022.27.6.2200042",     country = "Netherlands", subgroup = "Omicron" ),
  
  d36 = list(
    # Fig 2. Contact tracing, Japan, Jan 2022. L452R-confirmed Omicron BA.1, n=77.
    # Single known exposure date -> point observations, no interval censoring.
    # Digitised from density histogram; n=77 and mean=3.04 verified against reported 3.03.
    # Supersedes summary entry d36 (mean=3.03, sd=1.33, n=77).
    freq_value = c(1,  2,  3,  4, 5, 6, 7),
    freq_count = c(7, 21, 28, 10, 6, 4, 1),
    source   = "Tanaka (2022), doi: 10.3390/ijerph19106330",
    country  = "Japan",
    subgroup = "Omicron"
  ),
  
  d39 = list( mean = 3.5,    sd = 1.4,               n = 22,   source = "Liu (2022), doi: 10.1016/j.onehlt.2022.100425",                  country = "South Korea", subgroup = "Omicron" ),
  d43 = list( median = 3,    Q1 = 2,    Q3 = 4,      n = 36,   source = "Zeng (2023), doi: 10.3201/eid2904.220854",                       country = "Singapore",   subgroup = "Omicron" ),
  d44 = list( mean = 3.27,   sd = 1.05,              n = 500,  source = "Xiong (2023), doi: 10.3201/eid2902.221243",                      country = "China",       subgroup = "Omicron" ),
  d45 = list( mean = 4.6,    sd = 2.1,               n = 52,   source = "Wei (2023), doi: 10.1111/irv.13097",                             country = "China",       subgroup = "Omicron" ),
  
  # --- Delta ---
  d33 = list( mean = 4.4,    sd = 2.5,               n = 255,  source = "Backer (2022), doi: 10.2807/1560-7917.ES.2022.27.6.2200042",     country = "Netherlands", subgroup = "Delta" ),
  d37 = list( median = 3,    min = 1,    max = 13,   n = 171,  source = "McAleavey (2022), doi: 10.1016/j.puhe.2022.06.023",              country = "Ireland",     subgroup = "Delta" ),
  d38 = list( mean = 4.4,    sd = 1.9,               n = 47,   source = "Zhang (2021), doi: 10.46234/ccdcw2021.148",                      country = "China",       subgroup = "Delta" ),
  d40 = list( mean = 6.5,    sd = 3.7,               n = 64,   source = "Liu (2022), doi: 10.1016/j.onehlt.2022.100425",                  country = "South Korea", subgroup = "Delta" ),
  d42 = list( median = 4,    Q1 = 3,    Q3 = 7,      n = 42,   source = "Zeng (2023), doi: 10.3201/eid2904.220854",                       country = "Singapore",   subgroup = "Delta" ),
  d46 = list( mean = 5.3,    sd = 3.35,              n = 71,   source = "Luo (2023), doi: 10.46234/ccdcw2023.011",                        country = "China",       subgroup = "Delta" ),
  
  # --- Alpha ---
  d35 = list(
    # Fig 3. Contact tracing, Japan, Apr-May 2021. N501Y-confirmed Alpha, n=51.
    # Single known exposure date -> point observations, no interval censoring.
    # Digitised from density histogram; n=51, mean=4.92 vs reported 4.94.
    # Supersedes summary entry d35 (mean=4.94, sd=2.19, n=51).
    freq_value = c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12),
    freq_count = c( 2, 5, 9, 9, 8, 6, 5, 3, 2,  1,  1),
    source   = "Tanaka (2022), doi: 10.3390/ijerph19106330",
    country  = "Japan",
    subgroup = "Alpha"
  ),
  
  # --- Beta ---
  d41 = list( median = 4.5,  min = 2,    max = 7,     n = 10,   source = "Investigation team (2021), doi: 10.2807/1560-7917.ES.2021.26.13.2100333",  country = "France", subgroup = "Beta" ),
  
  # --- Frequency table data ---
  d_bohmer_2020 = list(         # Table. N=12 cases, Bavaria Germany, Wildtype. Individual onset and infection date ranges extracted. 4 cases excluded: P6/P9/P11 (unknown exposure date), P15 (asymptomatic). Supersedes summary-stat d18 (same study).
    freq_lower = c(2, 3, 1, 2, 2, 4, 4, 2, 5, 4, 3, 7),
    freq_upper = c(3, 5, 1, 4, 2, 4, 4, 2, 5, 7, 6, 7),
    freq_count = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    source = "Böhmer (2020), doi: 10.1016/S1473-3099(20)30314-5",
    country = "Germany",
    subgroup = "Wildtype"
  )
)


# Dengue ---------------------------------------------------------------------

# Underlying data based on VBD Systematic Review by Rudolph et al (2014)
# published in AJTMH: https://doi.org/10.4269/ajtmh.13-0403

datasets_Dengue <- list(
  d1 = list(       #Experimental human inoculation study, Philippines. Cases 1-6 and 9-10: intravenous blood inoculation; Case 11: mosquito bite (approximate). Case 6 prolonged incubation attributed by authors to probable relative immunity.
    freq_value = c(3.750, 2.792, 2.750, 2.500, 4.167, 7.000, 3.458, 2.500, 3.667),
    freq_count = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
    country = "Philippines",
    subgroup = "inoculation",
    source = "Ashburn & Craig 1907, https://www.jstor.org/stable/30073165"
  ),
  d2 = list(            #Experimental human inoculation study, Sydney, Australia. Subcutaneous injection route throughout. Case 6 (15 days) excluded: congenital syphilis, authors noted very atypical incubation. Case 17 excluded: indeterminate onset.
    freq_value = c(8, 6.5, 10, 9, 7, 8, 5, 7, 4.792),
    freq_count = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
    country = "Australia",
    subgroup = "inoculation",
    source = "Cleland (1919), doi: 10.1017/S0022172400007476"
  ),
  d3 = list(            #Experimental human inoculation and mosquito transmission study, Sydney, Australia. Subcutaneous injection route for inoculation cases. Mosquito cases via Stegomyia fasciata bites. Cases 6 and 7 attributed to first injection by authors; Case 9 to second injection. Case 17 onset taken as fever onset (~5.75 days); Case 29 onset taken as first symptoms (4d 21h). Mosquito Case IV (Wm.) interval-censored (bitten on two successive days).
    freq_lower = c(5.833, 5.75, 4.875, 6.667, 6.583, 6.583, 7.875, 8.125, 7.833, 8.542, 9, 8.5,
                   8.208, 9.417, 7.75, 5.208),
    freq_upper = c(5.833, 5.75, 4.875, 6.667, 6.583, 6.583, 7.875, 8.125, 7.833, 8.542, 9, 8.5,
                   8.208, 9.417, 7.75, 6.375),
    freq_count = rep(1, 16),
    country = "Australia",
    subgroup = "inoculation",
    source = "Cleland (1918), doi: 10.1017/S0022172400006690"
  ),
  d4 = list(
    freq_lower   = c( 4,  5,  6,  7,  8,  9, 10),
    freq_upper   = c( 5,  6,  7,  8,  9, 10, 11),
    freq_count   = c(11, 18,  7,  6,  3,  1,  1),
    country      = "Philippines",
    subgroup     = "mosquito bite",
    source       = "Siler (1926), https://archive.org/details/act3868.0029.001.umich.edu"
  ),
  d5 = list(
    freq_value   = c(6.50, 6.75, 7.25),
    freq_count   = c(1,   1,   1),
    country      = "Philippines",
    subgroup     = "inoculation",
    source       = "Siler (1926), https://archive.org/details/act3868.0029.001.umich.edu"
  ),
  d6 = list(       # exposure via single night under netting with infected Culex fatigans mosquitoes; exact bite time within night unknown
    freq_lower = c(4, 5, 6),
    freq_upper = c(4, 5, 6),
    freq_count = c(1, 1, 1),
    country    = "Lebanon",
    subgroup   = "inoculation",
    source     = "Graham (1903), https://wellcomecollection.org/works/psrzwk8h"
  ),
  d7 = list(       # single-bite cases (Walls, Vorhies, Herron) have equal bounds; re-bitten cases (Swartz, Mackinson, Werner, Colegrove, Bibeault) interval-censored with lower = days from second bite to onset, upper = days from first bite to onset
    freq_lower = c(5,  5,  5,  6,  7,  7,  9,  9),
    freq_upper = c(5,  15, 16, 6,  7,  17, 9,  20),
    freq_count = c(1, 1, 1, 1, 1, 1, 1, 1),
    country    = "Philippines",
    subgroup   = "mosquito bite",
    source     = "Schule PA (1928), doi: 10.4269/ajtmh.1928.s1-8.203"
  ),
  d8 = list(       # 1922 Texas epidemic; single bite event per case; IP measured to prodromal onset
    freq_value = c(4.08, 5.75, 6.17, 6.50),
    freq_count = c(1, 1, 1, 1),
    country    = "USA",
    subgroup   = "mosquito bite",
    source     = "Chandler (1923), doi: 10.4269/ajtmh.1923.s1-3.233"
  ),
  d9 = list(       # 1922 Texas epidemic; IP measured to prodromal onset; different route from mosquito bite cases
    freq_value = c(5.13, 5.79),
    freq_count = c(1, 1),
    country    = "USA",
    subgroup   = "inoculation",
    source     = "Chandler (1923), doi: 10.4269/ajtmh.1923.s1-3.233"
  )
)

# Yellow Fever ---------------------------------------------------------------------

# Underlying data based on VBD Systematic Review by Rudolph et al (2014)
# published in AJTMH: https://doi.org/10.4269/ajtmh.13-0403

datasets_YFV <- list(
  # the below is also contained in https://archive.org/details/yellowfeverepide00cart/page/6/mode/2up
  # there is data on 12 further cases which can't currently be accessed.
  # Reed (1902), Table II — Camp Lazear, mosquito bite
  # Incubation periods given to the half-hour; encoded as exact values
  d1 = list(
    freq_value = c(3.3958, 5.7083, 3.4792, 3.8125, 3.9583,
                   3.9375, 3.9792, 3.1042, 3.2500, 2.9167),
    freq_count = rep(1, 10),
    country    = "Cuba",
    subgroup   = "mosquito bite",
    source     = "Reed (1902), PMID: 20474139"
  ),
  # Reed (1902), Table III — Guitéras inoculation station, Havana
  # Incubation periods given as days + hours
  d2 = list(
    freq_value = c(3.4167, 4.2083, 3.1250, 5.1250,
                   3.7917, 3.8750, 5.8750, 3.0000),
    freq_count = rep(1, 8),
    country    = "Cuba",
    subgroup   = "mosquito bite",
    source     = "Reed (1902), PMID: 20474139"
  ),
  # Reed (1902), Table I — blood injection experiments
  # Only attack dates recorded (integer days); encoded as unit-width windows
  # Cases IV (2d) and XII (1d) flagged as uncertain due to implausibly short intervals
  # subcutaneous blood injection (0.5-2cc venous blood); integer-day resolution only, encoded as unit-width windows; cases IV (2d) and XII (1d) implausibly short, possibly date-recording artefacts
  d3 = list(
    freq_lower = c(3.5, 2.5, 1.5, 2.5, 4.5, 3.5, 3.5, 0.5),
    freq_upper = c(4.5, 3.5, 2.5, 3.5, 5.5, 4.5, 4.5, 1.5),
    freq_count = rep(1, 8),
    country    = "Cuba",
    subgroup   = "inoculation",
    source     = "Reed (1902), PMID: 20474139"
  )
)


# Flu ---------------------------------------------------------------------

# USE SUBGROUP FOR WHICH FLU


# Anthrax ---------------------------------------------------------------------

# Note the importance of dose-response model (can cite https://doi.org/10.1073/pnas.0509551103) 
