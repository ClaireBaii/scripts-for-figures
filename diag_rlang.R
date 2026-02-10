cat("LibPaths:\n")
print(.libPaths())
tryCatch({
    cat("rlang path: ", find.package("rlang"), "\n")
    cat("rlang version: ", as.character(packageVersion("rlang")), "\n")
}, error = function(e) cat("rlang not found or error:", e$message, "\n"))
