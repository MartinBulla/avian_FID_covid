require(rmarkdown)

rmarkdown::render("R/REV_ms_output.R", output_dir = "Outputs", output_file = "REV_ms_output_v3.html")
rmarkdown::render("R/REV_analysis.R", output_dir = "Outputs", output_file = "REV_analysis_v2.html")
rmarkdown::render("R/MET_google-mob_out.R", output_dir = "Outputs", output_file = "EXP_GoogleMobilitty.html")
