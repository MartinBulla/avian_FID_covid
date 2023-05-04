
# number of obseervers per country
table(d$Country, d$Observer)

# number of observers per site in Australia
da <- d[Country %in% "Australia"]
table(da$Observer, da$IDLocality)

# number of sites
length(unique(d$IDLocality))  

# number of Ausstralian sites
length(unique(da$IDLocality))

# number of unique day-years and hour-day-years
length(unique(paste(d$Day, d$Year))) 

db <- d[Covid == 0]
length(unique(paste(db$Day, db$Year)))
length(unique(paste(db$Day, db$Year, db$Hour)))

dc <- d[Covid == 1]
length(unique(paste(dc$Day, dc$Year)))
length(unique(paste(dc$Day, dc$Year, dc$Hour)))

# species differences for before and during covid periods
spb <- d[Covid == 0, unique(Species)]
spd <- d[Covid == 1, unique(Species)]
spd[!spd %in% spb]