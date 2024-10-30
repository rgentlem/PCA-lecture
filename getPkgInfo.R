## from ChatGPT - a helper script to write out all packages being used

# Get package info from session
sI <- sessionInfo()

pkg_info = sI$otherPkgs
# Extract package names and versions
pkg_df <- data.frame(
  Package = names(pkg_info),
  Version = sapply(pkg_info, function(pkg) pkg$Version)
)

# Save to a CSV file
write.csv(pkg_df, "package_versions.csv", row.names = FALSE)
