suppressMessages(library(cbioportalR))
set_cbioportal_db("public")
samps <- available_samples("brca_aurora_2023")
cat("n samples:", nrow(samps), "\n")
print(head(samps$sampleId, 5))
sid <- samps$sampleId[1]
m <- try(get_mutations_by_sample("brca_aurora_2023", sid), silent = TRUE)
if (inherits(m, "try-error")) {
  cat("by-sample failed:", as.character(m), "\n")
} else {
  cat("rows:", nrow(m), "| cols:", paste(names(m), collapse = ","), "\n")
}
