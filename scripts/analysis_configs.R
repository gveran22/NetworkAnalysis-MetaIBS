analysis_configs <- list(
  list(agg_level = "Order", filtTax = "numbSamp", filtTaxPar = list(numbSamp = 5), filtSamp = "totalReads", filtSampPar = list(totalReads = 500)),
  list(agg_level = "Family", filtTax = "numbSamp", filtTaxPar = list(numbSamp = 5), filtSamp = "totalReads", filtSampPar = list(totalReads = 500)),
  list(agg_level = "Genus", filtTax = "highestVar", filtTaxPar = list(highestVar = 100), filtSamp = "totalReads", filtSampPar = list(totalReads = 500))
)