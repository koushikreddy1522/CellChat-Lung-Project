
# Download the CellChat human lung dataset
# Source: Mendeley Data
# DOI: 10.17632/7r2cwbw44m.2

# Install jsonlite if needed
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

# Create the data folder
dir.create("datasets", showWarnings = FALSE)

# Dataset information
dataset_id <- "7r2cwbw44m"
version <- 2

# Get the list of dataset files from Mendeley
api_url <- paste0(
  "https://api.data.mendeley.com/datasets/public/",
  dataset_id, "/files?version=", version
)

files <- jsonlite::fromJSON(
  api_url,
  simplifyVector = FALSE
)

# Download the two required files separately
required_files <- c(
  "Counts_lung_cells.csv",
  "Metadata_lung_cells_update.csv"
)

for (filename in required_files) {
  
  file_info <- Filter(
    function(x) identical(x$filename, filename),
    files
  )
  
  if (length(file_info) == 0) {
    stop(paste("File not found:", filename))
  }
  
  file_id <- file_info[[1]]$id
  
  details_url <- paste0(
    "https://api.data.mendeley.com/datasets/public/",
    dataset_id, "/files/", file_id,
    "?version=", version
  )
  
  details <- jsonlite::fromJSON(details_url)
  
  download_url <- details$content_details$download_url
  
  if (is.null(download_url)) {
    stop(paste("Download URL unavailable for:", filename))
  }
  
  download.file(
    download_url,
    file.path("datasets", filename),
    mode = "wb"
  )
  
  cat("Downloaded:", filename, "\n")
}
