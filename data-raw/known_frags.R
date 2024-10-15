## code to prepare `known_frags` dataset goes here
known_frags <- unlist(list("M+H"=0,"M+"=-1,"M-H"=-2,"M+H2O-CH4"=+2))
usethis::use_data(known_frags, overwrite = TRUE)
