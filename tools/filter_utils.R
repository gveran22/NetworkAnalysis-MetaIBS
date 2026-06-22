
# **********************************
# Purpose: Help Functions- Filtering
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# **********************************

resolve_filter_params <- function(physeq_obj, filtTax, filtTaxPar, filtSamp, filtSampPar) {
  
  n_samples_ds <- nsamples(physeq_obj)
  samp_reads <- sample_sums(physeq_obj)
  
  filtTaxPar_resolved <- filtTaxPar
  filtSampPar_resolved <- filtSampPar
  
  # Resolver numbSamp proporcional
  if (!is.null(filtTax) && !"none" %in% filtTax && "numbSamp" %in% filtTax) {
    if (!is.null(filtTaxPar) &&!is.null(filtTaxPar$numbSamp_prop)) {
      filtTaxPar_resolved <- list(
        numbSamp = ceiling(filtTaxPar$numbSamp_prop * n_samples_ds)
      )
    }
  }
  
  # Resolver totalReads como percentil de profundidad por muestra
  if (!is.null(filtSamp) && !"none" %in% filtSamp && "totalReads" %in% filtSamp) {
    if (!is.null(filtSampPar) && !is.null(filtSampPar$totalReads_quant)) {
      q <- filtSampPar$totalReads_quant
      
      if (q == 0) {
        thr <- 0
      } else {
        thr <- unname(stats::quantile(samp_reads, probs = q, na.rm = TRUE, type = 7))
        thr <- floor(thr)
      }
      
      filtSampPar_resolved <- list(
        totalReads = thr
      )
    }
  }
  
  list(
    filtTaxPar_resolved = filtTaxPar_resolved,
    filtSampPar_resolved = filtSampPar_resolved
  )
}

# Filter samples by totalReads
filter_samples_totalReads <- function(physeq, totalReads_quant) {
  
  samp_reads <- phyloseq::sample_sums(physeq)
  
  if (totalReads_quant == 0) {
    thr <- 0
  } else {
    thr <- unname(stats::quantile(samp_reads, probs = totalReads_quant,
                                  na.rm = TRUE, type = 7))
    thr <- floor(thr)
  }
  
  keep <- samp_reads >= thr
  removed_samples <- names(samp_reads)[!keep]
  
  physeq_filt <- phyloseq::prune_samples(keep, physeq)
  
  list(
    physeq = physeq_filt,
    removed = removed_samples,
    kept = names(samp_reads)[keep],
    totalReads_used = thr
  )
}


# Filter taxa by numbSamp
filter_taxa_numbSamp <- function(physeq, numbSamp_prop) {
  
  n_samples_ds <- nsamples(physeq)
  otu <- phyloseq::otu_table(physeq)
  
  if (phyloseq::taxa_are_rows(physeq)) {
    prev <- rowSums(otu > 0)
  } else {
    prev <- colSums(otu > 0)
  }
  
  numbSamp = ceiling(numbSamp_prop * n_samples_ds)
  
  keep <- prev >= numbSamp
  
  removed_taxa <- names(prev)[!keep]
  
  physeq_filt <- phyloseq::prune_taxa(keep, physeq)
  
  list(
    physeq = physeq_filt,
    removed = removed_taxa,
    kept = names(prev)[keep],
    numbSamp_used = numbSamp
  )
}