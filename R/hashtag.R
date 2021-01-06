#' Compute log Fold Difference between highest hash tag and Closest Contender
#' (fdcc)
#'
#'
compute_fdcc <- function(hash_tag_counts) {
  if (is.null(hash_tag_counts) || length(unique(hash_tag_counts)) == 2) {
    return(NA)
  }
  count_order <- order(hash_tag_counts)
  max_entry <- hash_tag_counts[last(count_order)]
  sec_best <- hash_tag_counts[count_order[length(hash_tag_counts) - 1]]
  return(log2(max_entry + 1) - log2(sec_best + 1))
}


#' Compute Pielou evennness of distribution encoded in a vector
#'
#'
compute_evenness <- function(x, ignore_zeros = F) {
  stopifnot(class(x) == 'numeric' || class(x) == 'integer')
  if (ignore_zeros) {
    x <- x[x != 0]
  }
  x_norm <- x / sum(x, na.rm = T)
  ## Divide observed entropy by maximally obtainable entropy for the number of
  ## observed classes
  -sum(log((x_norm)^(x_norm))) / log(length(x))
}


extract_hashtags_from_cellranger <- function(
  data_dir = '/DATA/users/m.slagter/MirjamHoekstra/raw_exp_5310',
  gene_column = 2) {
  data_obj <- Seurat::Read10X(data.dir = data_dir, gene.column = gene_column,
                              strip.suffix = T)
  ht_i <- length(data_obj)
  hashtag_counts <- data_obj[[ht_i]]
  return(hashtag_counts)
}


compute_hashtag_stats <- function(hashtag_counts, fd_thresh = 2, 
                                  evenness_thresh = .5, read_thresh = 100) {
  stats <-
    data.frame(
      fd_cc = apply(hashtag_counts, 2, compute_fdcc),
      hashtag_evenness = apply(hashtag_counts, 2, compute_evenness),
      dominant_hashtag = apply(hashtag_counts, 2, which.max),
      total_hashtag_reads = apply(hashtag_counts, 2, sum, na.rm = T)
    ) %>%
    cbind(t(as.data.frame(hashtag_counts))) %>%
    tibble::rownames_to_column('sample_id') %>%
    dplyr::mutate(fd_crit = is.na(fd_cc) | fd_cc >= fd_thresh) %>%
    dplyr::mutate(evenness_crit = hashtag_evenness <= evenness_thresh) %>%
    dplyr::mutate(read_crit = total_hashtag_reads >= read_thresh)
  return(stats)
}


filter_seurat_stats <- function(seurat_object, stats) {
  allowed_sample_ids <- stats %>% 
    dplyr::filter(fd_crit == T & evenness_crit == T & read_crit == T) %>%
    pull(sample_id) %>%
    intersect(colnames(seurat_object))
  ret_val <- tryCatch(seurat_object[, allowed_sample_ids], 
                      error = function(e) { print(e); NULL }) 
  return(ret_val)
}
