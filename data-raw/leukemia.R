library(R.matlab)

# Load the MATLAB data file
leukemia <- readMat("data-raw/leukemia-raw.mat")
names(leukemia) <- c("x", "y")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(leukemia, overwrite = TRUE, compress = "xz")
