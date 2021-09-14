umf.df <- as(umf, "data.frame")

rugr.y <- ifelse(cbind(grouse1=apply(umf.df[,1:10]>0, 1, any, na.rm=TRUE),
                       grouse2=apply(umf.df[,11:20]>0, 1, any, na.rm=TRUE),
                       grouse3=apply(umf.df[,21:30]>0, 1, any, na.rm=TRUE)), 1L, 0L)

rugr.out <- with(umf.df,
                 data.frame(Route, Year, PointID, rugr.y, elev, fc, doy.1, doy.2, doy.3))


rugr <- rugr.out

write.csv(rugr, file="grouse_data.csv", row.names=TRUE)
