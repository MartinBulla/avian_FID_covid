require(rmarkdown)

rmarkdown::render("R/REV_analysis.R", output_dir = "Outputs", output_file = "REV_analysis.html")
rmarkdown::render("R/MET_google-mob_out.R", output_dir = "Outputs", output_file = "EXP_GoogleMobilitty.html")
