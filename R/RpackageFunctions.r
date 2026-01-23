


#' Convert SIRC results to DNAStringSet
#'
#' @description
#' Extracts SIRC (Short Inverted Repeat Cassette) sequences from a genome FASTA file
#' and converts them into a \code{\link[Biostrings]{DNAStringSet}} object.
#'
#' @param results R list of results obtained via \code{\link{SIRCFinder}} function.
#' @param selection Type of sequences to extract:
#'   \itemize{
#'     \item \code{"Cassettes"} - full SIRC sequences (default)
#'     \item \code{"DR"} - direct repeat elements only
#'     \item \code{"Spacers"} - spacer elements only
#'   }
#'   When selecting \code{"DR"} or \code{"Spacers"}, sequence IDs will include
#'   additional suffix \code{";N"} where N is the element number within the cassette.
#' @param path_to_genome Path to the genome FASTA file.
#'
#' @details
#' This function reads the genome FASTA file and extracts SIRC sequences based on
#' genomic coordinates provided in the \code{results} object.
#'
#' @return
#' A \code{\link[Biostrings]{DNAStringSet}} object containing the requested sequences.
#'
#' @section Sequence ID format:
#' \itemize{
#'   \item For \code{selection = "Cassettes"}: \code{[chromosome]_[cassette_id]}
#'   \item For \code{selection = "DR"} or \code{"Spacers"}: \code{[chromosome]_[cassette_id];[element_number]}
#' }
#'
#' @author
#' Igor V Gorbenko \email{gorbenko@@sifibr.irk.ru}
#' Karolina M. Zverintseva \email{zverincewa@@mail.ru}
#'
#' @references
#' Reference to be added after manuscript submission to bioRxiv.
#'
#' @examples
#' \dontrun{
#' # Assuming you have SIRC results and a genome file
#' sirc_results <- SIRCFinder(...)
#'
#' # Extract full cassettes
#' cassettes <- convert_to_DNA(sirc_results,
#'                             selection = "Cassettes",
#'                             path_to_genome = "genome.fasta")
#'
#' # Extract only direct repeats
#' drs <- convert_to_DNA(sirc_results,
#'                       selection = "DR",
#'                       path_to_genome = "genome.fasta")
#' }
#'
#' @export
convert_to_DNA <- function(results, selection = "Cassettes", path_to_genome){
  require(msa, quietly = T)
  require(data.table, quietly = T)


  make_all_unique_dt <- function(x, sep = ";") {
    dt <- data.table(x = x, idx = seq_along(x))
    dt[, count := seq_len(.N), by = x]
    dt[, result := paste0(x, sep, count)]
    dt[order(idx)]$result
  }


  gr <- convert_to_granges(res)

  if(selection == "Cassettes"){
    gr <- gr %>% filter(type == "cassette")
  } else if(selection == "DR"){
    gr <- gr %>% filter(type == "DR")
  } else if(selection == "Spacers") {
    gr <- gr %>% filter(type == "SP")
  } else {return(print("selection must be one of: Cassettes, DR, Spacers !"))}

  gr2 <- as.data.frame(gr)
  nms <- paste0(gr2$seqnames, "_", gr2$cassette_id)

  if(selection != "Cassettes"){
    nms <- make_all_unique_dt(nms)}

  gen <- readDNAStringSet(path_to_genome)
  names(gen) <- stringr::str_split_fixed(names(gen), " ",2)[,1]
  out <- gen[gr]
  names(out) <- nms
  return(out)
}



#' Convert R list of SIRC results to GenomicRanges
#'
#' @description
#' Converts the R list of SIRC detection results into a
#' \code{\link[GenomicRanges]{GRanges}} object, either in full form
#' (containing all DR and spacer positions) or short form (only cassette
#' positions).
#'
#' @param rlist R list of results obtained via \code{\link{SIRCFinder}} function.
#' @param type Type of output: either \code{"full"} for output containing
#'   all SIRC elements positions (cassettes, DRs, and spacers), or
#'   \code{"short"} for output containing only full cassette positions
#'   (default: \code{"full"}).
#'
#' @details
#' This function extracts genomic positions from the SIRC results list and
#' converts them into a \code{GRanges} object. The conversion includes:
#' \itemize{
#'   \item Sequence names (chromosomes/contigs)
#'   \item Start and end positions for each element
#'   \item Element type: "cassette", "DR", or "SP" (spacer)
#'   \item Cassette ID linking elements to their parent cassette
#'   \item Additional metadata: number of spacers per cassette, DR length,
#'         and DR consensus sequence
#' }
#'
#' For \code{type = "full"}, the output includes all individual elements
#' (cassette boundaries, direct repeats, and spacers). For \code{type = "short"},
#' only the boundaries of complete cassettes are returned.
#'
#' @return
#' A \code{\link[GenomicRanges]{GRanges}} object containing genomic ranges
#' of SIRC elements with the following metadata columns:
#' \itemize{
#'   \item \code{type}: Type of element ("cassette", "DR", or "SP")
#'   \item \code{cassette_id}: Integer ID of the parent cassette
#'   \item \code{spacers_number}: Number of spacers in the cassette
#'   \item \code{dr_length}: Length of direct repeats in base pairs
#'   \item \code{DR_cons}: Consensus sequence of direct repeats
#' }
#'
#' @author
#' Igor V Gorbenko \email{gorbenko@@sifibr.irk.ru}
#' Karolina M. Zverintseva \email{...}
#'
#' @references
#' Gorbenko et al. (2026). SIRCFinder: detection of Short Inverted Repeat
#' Cassettes in eukaryotic genomes. bioRxiv. (Reference to be added after
#' manuscript submission)
#'
#' @examples
#' \dontrun{
#' # Get SIRC results
#' results <- SIRCFinder("genome.fasta")
#'
#' # Convert to full GRanges with all elements
#' gr_full <- convert_to_granges(results, type = "full")
#'
#' # View the structure
#' gr_full
#'
#' # Get only cassette boundaries
#' gr_short <- convert_to_granges(results, type = "short")
#'
#' # Count elements by type in full version
#' table(gr_full$type)
#' }
#'
#' @export
convert_to_granges <- function(rlist, type = "full") {
  require(plyranges)
  convert_internal <- function(rlist){total_elements <- sum(sapply(rlist, function(x) {
    n <- length(x$positions)
    1 + n + max(0, n - 1)
  }))


  seqnames <- character(total_elements)
  starts <- integer(total_elements)
  ends <- integer(total_elements)
  types <- character(total_elements)
  cassette_ids <- integer(total_elements)

  idx <- 1
  for (i in seq_along(rlist)) {
    cassette <- rlist[[i]]
    seqname <- cassette$seqname
    positions <- sort(cassette$positions)
    dr_length <- cassette$dr_length
    n <- length(positions)

    seqnames[idx] <- seqname
    starts[idx] <- positions[1]
    ends[idx] <- positions[n] + dr_length - 1
    types[idx] <- "cassette"
    cassette_ids[idx] <- i
    idx <- idx + 1

    for (j in seq_len(n)) {
      seqnames[idx] <- seqname
      starts[idx] <- positions[j]
      ends[idx] <- positions[j] + dr_length - 1
      types[idx] <- "DR"
      cassette_ids[idx] <- i
      idx <- idx + 1
    }

    if (n > 1) {
      for (k in 1:(n - 1)) {
        spacer_start <- positions[k] + dr_length
        spacer_end <- positions[k + 1] - 1

        if (spacer_start <= spacer_end) {
          seqnames[idx] <- seqname
          starts[idx] <- spacer_start
          ends[idx] <- spacer_end
          types[idx] <- "SP"
          cassette_ids[idx] <- i
          idx <- idx + 1
        }
      }
    }
  }

  if (idx <= total_elements) {
    seqnames <- seqnames[1:(idx-1)]
    starts <- starts[1:(idx-1)]
    ends <- ends[1:(idx-1)]
    types <- types[1:(idx-1)]
    cassette_ids <- cassette_ids[1:(idx-1)]
  }

  as_granges(data.frame(
    seqnames = seqnames,
    start = starts,
    end = ends,
    type = types,
    cassette_id = cassette_ids,
    stringsAsFactors = FALSE
  ))


  }
  temp <- convert_internal(rlist)
  idx <- temp$cassette_id
  spacers_numbers <- sapply(rlist, function(x) length(x$positions))
  dr_lengths <- sapply(rlist, function(x) x$dr_length)
  DR_conss <- sapply(rlist, function(x) x$consensus)

  temp <- temp %>% mutate(spacers_number = spacers_numbers[idx],
                          dr_length = dr_lengths[idx],
                          DR_cons = DR_conss[idx])

  temp <- as.data.frame(temp)
  temp <- temp[order(temp$cassette_id, temp$start),]
  temp <- as_granges(temp)

  if(type == "short"){
    temp <- temp %>% filter(type == "cassette")
  }

  return(temp)

}


#' Remove scaffolds from assembly fasta files
#'
#' @description
#' Processes multifasta files by removing sequences whose names contain
#' the word "scaffold". Useful for filtering out scaffold sequences from
#' genome assemblies before SIRC analysis.
#'
#' @param files Character vector of paths to FASTA files for processing.
#'   Each file will be processed independently.
#'
#' @details
#' The function reads each input FASTA file, identifies sequences with
#' names containing the substring "scaffold" (case-insensitive), removes
#' these sequences, and writes the filtered sequences to new files.
#' Output files are named by appending ".processed.fa" to the original
#' filenames.
#'
#' @return
#' The function does not return a value but creates processed files on disk.
#' For each input file, a corresponding output file with extension
#' ".processed.fa" is created in the same directory.
#'
#' @note
#' The function performs case-insensitive matching for "scaffold".
#' Sequences with names like "scaffold_001", "Scaffold", "SCAFFOLD42", etc.
#' will all be removed.
#'
#' @author
#' Igor V Gorbenko \email{gorbenko@@sifibr.irk.ru}
#' Karolina M. Zverintseva \email{...}
#'
#' @references
#' Gorbenko et al. (2026). SIRCFinder: detection of Short Inverted Repeat
#' Cassettes in eukaryotic genomes. bioRxiv. (Reference to be added after
#' manuscript submission)
#'
#' @examples
#' \dontrun{
#' # Process single file
#' unscaffold("genome_assembly.fasta")
#' # Output: genome_assembly.fasta.processed.fa
#'
#' # Process multiple files
#' unscaffold(c("assembly1.fa", "assembly2.fa", "assembly3.fa"))
#' # Output: assembly1.fa.processed.fa, assembly2.fa.processed.fa, etc.
#' }
#'
#' @export
unscaffold <- function(files){
  require(msa, quietly = T)
  for(i in seq_along(files)){
    cat(paste("cleaning", files[i], "..."), "\r")
    mul <- readDNAStringSet(files[i])
    mul <- mul[!grepl("scaffold", names(mul))]

    writeXStringSet(mul, paste0(files[i], ".processed.fa"))
  }
}
