library(odbc)
library(DBI)
file_path <- "C:/Users/Visitor/Desktop/Coweeta_Data/SouthernApps.accdb"
accdb_con <- dbConnect(drv = odbc(), 
                       .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",file_path,";"))


## ------------- Format the net data --------------------
library(ggplot2)
library(plyr)
neteffort <- dbReadTable(accdb_con , "tbl_NetEffort_day")
neteffort$Year <- format(neteffort$NettingDate, format = "%Y")

#neteffort <- neteffort[neteffort$Year == "2018",]



neteffort_day <- dbReadTable(accdb_con , "tbl_NetEffort_IndNet")
neteffort_day <- neteffort_day[neteffort_day$Array.Type == "Grid",]
neteffort_day$o <- 1
neteffort_day$p <- 0
neteffort_day[!is.na(neteffort_day$AddPlayback),"p"] <- 1
nets <- dbReadTable(accdb_con , "tbl_Nets")  
head(nets)
tail(neteffort_day)

library(tidyr)

nets_wide1 <- pivot_wider(data = neteffort_day, id_cols = c(NetDate, Net), 
                         names_from = NetDate, values_from = o,
                         values_fn = mean, names_prefix = "open")
nets_wide2 <- pivot_wider(data = neteffort_day, id_cols = c(NetDate, Net), 
                         names_from = NetDate, values_from = p,
                         values_fn = mean, names_prefix = "PB")
nets_wide <- cbind(nets_wide1, nets_wide2)
nets_wide <- nets_wide[,-c(190, 191, 2)]

nets_wide$site <- nets[match(nets_wide$Net, nets$NetID, nomatch = NA), "Site"]
nets_wide$UTM_W <- nets[match(nets_wide$Net, nets$NetID, nomatch = NA), "W_UTM"]
nets_wide$UTM_N <- nets[match(nets_wide$Net, nets$NetID, nomatch = NA), "N_UTM"]
nets_wide <- nets_wide[nets_wide$site %in% unique(nets_wide$site)[c(2:18, 20)],]
nets_wide[is.na(nets_wide)] <- 0

write.csv(nets_wide, "Constant_Effort_Nets.csv")
