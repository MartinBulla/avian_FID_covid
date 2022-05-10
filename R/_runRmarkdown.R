require(rmarkdown)

rmarkdown::render('R/vis_raw.R', output_dir = 'Reports',  output_file = 'FID_before-during_COVID_raw.html')