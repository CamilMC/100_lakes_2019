id_niva <- readxl::read_excel("7.1000_lakes/id_niva.xlsx")
id_cba <- readxl::read_excel("5.100_lakes/100lakesData.xlsx") %>% distinct(lake_id,.keep_all = T)
nve_number <- readxl::read_excel("5.100_lakes/100lakes_gases.xlsx") %>% select(c("Lake_ID","NVE_number")) %>% 
  distinct(Lake_ID,.keep_all = T)

id_key <- merge(id_cba,id_niva,by.x="lake_id",by.y="Lake_ID") %>% 
  merge(nve_number,by.x = "lake_id", by.y = "Lake_ID")
id_key$temp <- NULL
id_key$cond. <- NULL
id_key$pH_sampling <- NULL
id_key$week <- NULL
id_key$PROJECT <- NULL
id_key$TEXT_ID_1 <- NULL
id_key$TEXT_ID_2 <- NULL
id_key$LAKE_NAME <- NULL
id_key$NAME <- NULL
names(id_key) <- c("Lake_ID","Lake_name","Long","Lat","CBA_sample_date","NIVA_station_ID","NIVA_sample_date","NVE_number")
id_key <- id_key %>% select("Lake_ID","NIVA_station_ID","NVE_number","Lake_name","Long","Lat","CBA_sample_date","NIVA_sample_date")
writexl::write_xlsx(id_key,"ID_key.xlsx")
