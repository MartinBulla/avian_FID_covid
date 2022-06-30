---------------------------------------------------------------------------------------

Description of the primary data from: 

"Urban birds' flight responses were unaffected by the COVID-19 shutdowns"

by Peter Mikula, Martin Bulla, Daniel Blumstein, Yanina Benedetti, Kristina Floigl, 
Jukka Jokimäki, Marja-Liisa Kaisanlahti-Jokimäki, Gábor Markó, Federico Morelli, 
Anders Pape Moller, Anastasiia Siretckaia, Sára Szakony, Mike Weston, Farah Abou Zeid, 
Piotr Tryjanowski, Tomáš Albrecht

---------------------------------------------------------------------------------------

WHEN USING the data, PLEASE CITE both the original paper and this dataset:

Bulla, M. & Mikula, P. (2022). Supporting information for 'Urban birds' flight responses 
were unaffected by the COVID-19 shutdowns'. Open Science Framework 
https://doi.org/10.17605/OSF.IO/WUZH7.

---------------------------------------------------------------------------------------

1. data.txt
the primary tab-delimited table containing data on the avian tolerance towards humans 
(measured as the flight initiation distance) and a set of predictors in urban birds

IDObs - unique identification of each observation
Country - the name of the country where the escape distance was collected
IDLocality - unique identification of sampled site
Lat - latitude of the is (IDLocality) in decimal degrees
Lon - longitude of the is (IDLocality) in decimal degrees
Covid - period  of the escape distance collection (0 = before the COVID-19 shutdowns; 1 = during the COVID-19 shutdowns)
Year - year of escape distance collection
Day - day of escape distance collection in relation to the start of the breeding season (Europe: Day 1 = 1 April; Australia: Day 1 = 15 August)
Hour - the hour of the day within which the escape was collected 
Temp - the ambient air temperature (°C) at the site during escape collection 
StringencyIndex - the strength of governmental measures rescaled to values from 0 to 100 (0 = no restrictions; 100 = strictest restrictions) and collected from Our World in Data database (https://ourworldindata.org/covid-stringency-index), based on data originally published in Hale et al. (2021)
Species - the scientific (Latin) name of the species as used in the BirdTree project (http://birdtree.org/)
FlockSize - the flock size
SD - the starting distance of observer (m)
FID - the flight initiation distance (m)
BodyMass - mean species body mass (g) estimated as the mean of female and male values from EltonTraits 1.0 database (Wilman et al. 2014)
Observer - the name of the observer