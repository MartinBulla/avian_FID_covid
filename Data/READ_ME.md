### Description of the primary data from:

[Urban birds' flight responses were unaffected by the COVID-19 shutdowns](https://doi.org/10.1101/2022.07.15.500232) 

by Peter Mikula, [Martin Bulla](https://martinbulla.github.io), Daniel T. Blumstein, Yanina Benedetti, Kristina Floigl, Jukka Jokimäki, Marja-Liisa Kaisanlahti-Jokimäki, Gábor Markó, Federico Morelli, Anders Pape Møller, Anastasiia Siretckaia, Sára Szakony, Michael A. Weston, Farah Abou Zeid, Piotr Tryjanowski & Tomáš Albrecht  

https://doi.org/10.1101/2022.07.15.500232  

***  

**WHEN USING** the data, **PLEASE CITE** both the [paper](https://doi.org/10.1101/2022.07.15.500232) and this repository (Martin Bulla, Peter Mikula, Daniel T. Blumstein, Yanina Benedetti, Kristina Floigl, Jukka Jokimäki, Marja-Liisa Kaisanlahti-Jokimäki, Gábor Markó, Federico Morelli, Anders Pape Møller, Anastasiia Siretckaia, Sára Szakony, Michael A. Weston, Farah Abou Zeid, Piotr Tryjanowski & Tomáš Albrecht (2023), *Supporting information for 'Urban birds’ flight responses were largely unaffected by the COVID-19 shutdowns'*, GitHub, [https://martinbulla.github.io/avian_FID_covid/](https://martinbulla.github.io/avian_FID_covid/)).  

***

### 1. data.txt   
The primary tab-delimited table containing data on the avian tolerance towards humans 
(measured as the flight initiation distance) and a set of predictors in urban birds  

- **IDObs** - unique identification of each escape distance trial (observation)  
- **Observer** - the name of the observer conducting the escape distancee trial  
- **Country** - the name of the country where the escape distance was collected  
- **IDLocality** - unique identification of sampled site  
- **Lat** - latitude of the IDLocality in decimal degrees  
- **Lon** - longitude of the IDLocality in decimal degrees  
- **Species** - the scientific (Latin) name of the species as used in the [BirdTree project](http://birdtree.org/)  
- **BodyMass** - mean species body mass in grams estimated as the mean of female and male values from [EltonTraits 1.0 database (Wilman et al. 2014)](https://doi.org/10.1890/13-1917.1)  
- **date_** - date of the escape distance trial in a format yyyy-mm-dd  
- **Year** - year of the escape distance trial 
- **Day** - day of the escape distance trial in relation to the start of the breeding season (Europe: Day 1 = 1 April; Australia: Day 1 = 15 August)  
- **Hour** - the hour of the day f the escape distance triall  
- **Temp** - the ambient air temperature in °C at the site during the escape distancee trial   
- **Covid** - was the escape distance trial before (0) or during the period of  the COVID-19 shutdowns (1)
- **StringencyIndex** - the strength of governmental measures rescaled to values from 0 to 100 (0 = no restrictions; 100 = strictest restrictions) and collected from [Our World in Data database](https://ourworldindata.org/covid-stringency-index), based on data originally published in [Hale et al. (2021)](https://doi.org/10.1038/s41562-021-01079-8)  
- **Human** - # of humans within 50 meter radius during the escape distance trial  
- **Height** - perch height of the bird in meteres during the escape distance trial  
- **FlockSize** - the flock size during eescape distance trial  
- **SD** - the starting distance of the observer in meters to the bird  
- **FID** - the flight initiation distance (escape distance) of the bird from the approaching observer in meters  

### 2. google_mobility.txt
The tab-delimited table containing [Google Mobility Reports](https://www.google.com/covid19/mobility/) (see there for data descriptoin); the variable of interest is "parks_percent_change_from_baseline"   

### 3. phylopic.txt   
The tab-delimited table containing names and IDs of bird silhouettes in [PhyloPic database](http://phylopic.org/)
- **Name** - the scientific (Latin) name of the species
- **Code** - Unique ID of each silhouette in PhyloPic database

### 4. taxonomy.txt
The tab-delimited table containing species names and taxonomic ranking of species used
in this study  
- **Species** -  the scientific (Latin) name of the species as used in the [BirdTree project](http://birdtree.org/)  
- **Family** - family affiliation of the species  
- **Order** - order affiliation of th species  

### 5. trees.tre
Text file containing 100 randomly sampled phylogenetic trees (Hackett backbone) from the [Birdtree project](http://birdtree.org/)