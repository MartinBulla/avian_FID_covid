require(rmarkdown)

rmarkdown::render("R/REV_ms_output.R", output_dir = "Outputs", output_file = "REV_ms_output_v12.html")

# END