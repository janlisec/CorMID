## code to prepare `precalc_idx` dataset goes here
precalc_idx <- lapply(1:5, function(n) { lapply(1:c(80,80,20,10,10)[n], function(N) { get_idx_mat(n, N) }) })
usethis::use_data(precalc_idx, overwrite = TRUE)
