

# Function to download all FASTQ files for a given BioProject ID (by ChatGPT 40)
#' Title
#'
#' @param bioproject_id
#' @param output_dir
#'
#' @return
#' @export
#'
#' @examples
#'


download_fastq_for_bioproject <- function(bioproject_id, output_dir = ".") {
  # Check if output directory exists, if not, create it
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Step 1: Get the list of SRA run IDs for the BioProject using EDirect
  message("Fetching SRA run IDs for BioProject: ", bioproject_id)
  edirect_cmd <- paste0("esearch -db sra -query '", bioproject_id, "' | efetch -format runinfo | cut -d',' -f1")
  sra_run_ids <- system(edirect_cmd, intern = TRUE)

  # Step 2: Filter the valid SRA run IDs (skip the header)
  sra_run_ids <- sra_run_ids[-1]  # Remove the first line (header)

  if (length(sra_run_ids) == 0) {
    stop("No SRA run IDs found for BioProject: ", bioproject_id)
  }

  message(length(sra_run_ids), " SRA runs found.")

  # Step 3: Download and convert each SRA run to FASTQ
  for (sra_id in sra_run_ids) {
    message("Processing SRA run: ", sra_id)

    # Prefetch the SRA file
    prefetch_cmd <- paste("prefetch", sra_id, "-O", output_dir)
    message("Running prefetch: ", prefetch_cmd)
    system2(command = "prefetch", args = c(sra_id, "-O", output_dir))

    # Convert the SRA file to FASTQ format using fastq-dump
    sra_file <- file.path(output_dir, paste0(sra_id, ".sra"))
    fastq_dump_cmd <- paste("fastq-dump", "--split-files", "--gzip", "-O", output_dir, sra_file)
    message("Running fastq-dump: ", fastq_dump_cmd)
    system2(command = "fastq-dump", args = c("--split-files", "--gzip", "-O", output_dir, sra_file))
  }

  message("All runs for BioProject ", bioproject_id, " have been downloaded and converted.")
}


# Function to fetch BioSamples table for a given BioProject using NCBI EDirect (by ChatGPT 40)
#' Title
#'
#' @param bioproject_id
#'
#' @return
#' @export
#'
#' @examples
#'


fetch_biosamples_for_bioproject <- function(bioproject_id) {

  # Step 1: Use esearch and efetch to retrieve BioSamples linked to the BioProject
  message("Fetching BioSamples for BioProject: ", bioproject_id)
  edirect_cmd <- paste0(
    "esearch -db biosample -query '", bioproject_id, "' | efetch -format docsum"
  )

  # Run the system command and capture the output
  biosample_output <- system(edirect_cmd, intern = TRUE)

  # Check if any results were returned
  if (length(biosample_output) == 0) {
    stop("No BioSamples found for BioProject: ", bioproject_id)
  }

  # Step 2: Parse the XML output from efetch
  require(XML)

  # Convert the XML output to a structured R list
  parsed_data <- xmlTreeParse(biosample_output, asText = TRUE, useInternalNodes = TRUE)
  root_node <- xmlRoot(parsed_data)

  # Step 3: Extract relevant information (e.g., BioSample accession, organism, title, etc.)
  biosample_list <- xpathApply(root_node, "//DocSum", function(node) {
    biosample_id <- xpathSApply(node, ".//Item[@Name='Accession']", xmlValue)
    organism <- xpathSApply(node, ".//Item[@Name='Organism']", xmlValue)
    title <- xpathSApply(node, ".//Item[@Name='Title']", xmlValue)

    # Create a list for each BioSample
    list(
      BioSample = biosample_id,
      Organism = organism,
      Title = title
    )
  })

  # Convert the list of BioSamples into a data frame
  biosample_df <- do.call(rbind, lapply(biosample_list, function(x) data.frame(x, stringsAsFactors = FALSE)))

  message(nrow(biosample_df), " BioSamples found for BioProject ", bioproject_id)

  return(biosample_df)
}

# Load necessary library
library(rentrez)
library(dplyr)
library(XML)

# Function to fetch metadata, including attributes, for a given BioProject
fetch_biosample_metadata <- function(bioproject_id) {

  # Step 1: Search for SRA records linked to the BioProject
  message("Searching SRA for BioProject: ", bioproject_id)
  search_results <- entrez_search(db = "sra", term = bioproject_id, retmax = 1000)

  if (length(search_results$ids) == 0) {
    stop("No SRA records found for BioProject: ", bioproject_id)
  }

  # Step 2: Fetch the summaries for each SRA record to get BioSample IDs
  message("Fetching SRA summaries to extract BioSample IDs...")
  sra_summaries <- entrez_summary(db = "sra", id = search_results$ids)
  biosample_ids <- sapply(sra_summaries, function(sra) sra$biosample)

  biosample_ids <- biosample_ids[!is.na(biosample_ids)]  # Remove NULL values

  if (length(biosample_ids) == 0) {
    stop("No BioSample IDs found for the SRA records.")
  }

  message("Found ", length(biosample_ids), " unique BioSample IDs.")

  # Step 3: Fetch detailed metadata for each BioSample
  message("Fetching detailed BioSample metadata...")
  biosample_metadata <- lapply(biosample_ids, function(biosample_id) {
    biosample_xml <- entrez_fetch(db = "biosample", id = biosample_id, rettype = "xml")

    # Parse the XML to extract attributes and metadata
    biosample_parsed <- xmlParse(biosample_xml)
    biosample_root <- xmlRoot(biosample_parsed)

    # Extract core metadata
    accession <- xmlValue(biosample_root[["Ids"]][["Id"]])
    organism <- xmlValue(biosample_root[["Organism"]][["OrganismName"]])
    submission_date <- xmlValue(biosample_root[["Submission"]][["Date"]])
    title <- xmlValue(biosample_root[["Title"]])

    # Extract attributes (these are key-value pairs)
    attributes <- xpathSApply(biosample_root, "//Attribute", xmlValue)
    attribute_names <- xpathSApply(biosample_root, "//Attribute", xmlGetAttr, "attribute_name")

    # Combine attributes into a named list
    attr_list <- setNames(as.list(attributes), attribute_names)

    # Combine all metadata in a single list
    metadata_list <- c(
      list(Accession = accession, Organism = organism, Submission_Date = submission_date, Title = title),
      attr_list
    )

    return(metadata_list)
  })

  # Step 4: Convert the list of metadata to a data frame
  message("Converting metadata to a data frame...")
  biosample_df <- bind_rows(lapply(biosample_metadata, as.data.frame))

  return(biosample_df)
}

