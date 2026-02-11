# Install mixOmics
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Set CRAN mirror to avoid timeouts
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})

if (!requireNamespace("mixOmics", quietly = TRUE)) {
    message("Installing mixOmics...")
    BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
} else {
    message("mixOmics is already installed.")
}

if (!requireNamespace("circlize", quietly = TRUE)) {
    message("Installing circlize...")
    install.packages("circlize")
}
