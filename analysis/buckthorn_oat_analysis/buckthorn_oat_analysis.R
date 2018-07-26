#This script is for spatial analysis of buckthorn and oat data
#Bulk raw data for Rhamnus cathartica (european buckthorn) was obtained from the following source on 05/18/2018: EDDMapS. 2018. Early Detection & Distribution Mapping System. The University of Georgia - Center for Invasive Species and Ecosystem Health. Available online at http://www.eddmaps.org/; 
#The raw bulk data was first counted by observation date and county (thanks to Rebekah D. Wallace at the Center for Invasive Species & Ecosystem Health for the database query)
#The county data was then manually curated to include records from publication and historical records that didn't have an observation date in the database (since it wasn't a direct field observation)

#Oat production data was downloaded from the USDA-ARS Cereal Disease Lab website: https://www.ars.usda.gov/midwest-area/stpaul/cereal-disease-lab/docs/small-grain-losses-due-to-rust/small-grain-losses-due-to-rust/
#The production data was extracted from the small grain losses due to rust files for 1990 and 2015 and compiled into a single file for analysis

#Change to the appropriate working directory
setwd("")

#Load appropriate libraries
library(dplyr)
library(ggplot2)
library(scales) #To add commas to labels if need be
library(grid)

bt_data <- read.csv("buckthorn_FirstDateStateCountyCount.csv", header = TRUE)
isolate_data <- read.csv("isolate_locations.csv", header = TRUE)

bt_data$state <- c("alabama", "alaska", "arizona", "arkansas", "california", "colorado", "connecticut", "delaware", "florida", "georgia", "hawaii", "idaho", "illinois", "indiana", "iowa", "kansas", "kentucky", "louisiana", "maine", "maryland", "massachusetts", "michigan", "minnesota", "mississippi", "missouri", "montana", "nebraska", "nevada", "new hampshire", "new jersey", "new mexico", "new york", "north carolina", "north dakota", "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina", "south dakota", "tennessee", "texas", "utah", "vermont", "virginia", "washington", "west virginia", "wisconsin", "wyoming")

#Which isolates are classified as BT, NO BT, or MN BT
bt_data$host_status <- ifelse(bt_data$PercentCounties < 10, bt_data$host_status <- "NO BT", bt_data$host_status <- "BT")
isolate_data$host_status <- bt_data[match(sub("\\d+(\\w{2}).+", "\\1", isolate_data$IsolateID), bt_data$StateName), 6]
isolate_data$host_status[isolate_data$Location == "St Paul BT, MN"] <- "MN BT"

states_map <- map_data("state")

#Make plot
map <- ggplot() +
  geom_map(data = bt_data, aes(map_id = state, fill = PercentCounties), color = "black", map = states_map) +
  expand_limits(x = states_map$long, y = states_map$lat) +
  scale_fill_gradient(low="white", high="darkblue", name="Percent of Counties") +
  theme_void() +
  geom_point(data=isolate_data, aes(x=Lon, y=Lat, color=factor(Year)), position = position_jitter(width = 0.25, height = 0.25), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("2015" = "red", "1990" = "black"), name = "Year")  +
  theme(legend.position="top", legend.text=element_text(size=12)) +
  guides(fill = guide_legend(title.position = "top"), color = guide_legend(title.position = "top"))

#Save plot
pdf("bt_isolate_map_no_inset.pdf", height = 4, width = 5)
map
dev.off()

#Oat analysis
#Note that only two states had an apparent increase in oat production because there was no recorded production in 1990, but production shown in 2015, and these are not shown on the map. These states are Massachusetts and Virginia. 
oat_data <- read.csv("oat_production_1990_2015.csv", header = TRUE)

oat_data$StateName <- c("alabama", "alaska", "arizona", "arkansas", "california", "colorado", "connecticut", "delaware", "florida", "georgia", "hawaii", "idaho", "illinois", "indiana", "iowa", "kansas", "kentucky", "louisiana", "maine", "maryland", "massachusetts", "michigan", "minnesota", "mississippi", "missouri", "montana", "nebraska", "nevada", "new hampshire", "new jersey", "new mexico", "new york", "north carolina", "north dakota", "ohio", "oklahoma", "oregon", "pennsylvania", "rhode island", "south carolina", "south dakota", "tennessee", "texas", "utah", "vermont", "virginia", "washington", "west virginia", "wisconsin", "wyoming")

map <- ggplot() +
  geom_map(data=oat_data, aes(map_id = StateName, fill = production_decline), color = "black", map = states_map) +
  expand_limits(x = states_map$long, y = states_map$lat) +
  scale_fill_gradient(low="lightgreen", high="red", name="Production Decline (1,000 of bushels)", breaks=c(10000, 20000, 30000, 40000), labels=c("10,000", "20,000", "30,000", "40,000")) +
  theme_void() +
  theme(legend.position="top", legend.text=element_text(size=12)) +
  guides(fill = guide_legend(title.position = "top"))

#Save plot
pdf("oat_production_map.pdf", height = 4, width = 5)
map
dev.off()


#In a previous iteration I tried to have Ramsey county as an inset to show where the MNBT isolates come from but decided against it for clarity.
#county_df <- map_data("county")
#ramsey_countyMap <- county_df %>% filter(region == "minnesota") %>% filter(subregion == "ramsey")

#Below is the code to make the plot with the county inset and states outlined by host status
# #Make plot
# map <- ggplot() +
#   geom_map(data=bt_data, aes(map_id = state, fill = PercentCounties, color = host_status), map = states_map) +
#   expand_limits(x = states_map$long, y = states_map$lat) +
#   scale_fill_gradient(low="white", high="darkblue", name="Percent of Counties") +
#   theme_void() +
#   geom_point(data=isolate_data[isolate_data$Location != "St Paul BT, MN",], aes(x=Lon, y=Lat, color=factor(Year)), position = position_jitter(width = 0.25, height = 0.25), alpha = 0.7, size = .08) +
#   geom_polygon(data=ramsey_countyMap, aes(long, lat), fill = "#56B4E9") +
#   # geom_map(data=ramsey_countyMap, aes(x=long, y=lat, map_id = subregion), map = ramsey_countyMap) +
#   scale_color_manual(values = c("NO BT" = "#E69F00", "BT" = "#D55E00", "2015" = "red", "1990" = "black", "#56B4E9"), name = "Year and Alternate Host Status", guide = FALSE)
# 
# #Inset map of Ramsey county
# ramsey_map <- ggplot() +
#   geom_polygon(data=ramsey_countyMap, aes(long, lat), fill = "#56B4E9") +
#   geom_point(data=isolate_data[isolate_data$Location == "St Paul BT, MN",], aes(x=Lon, y=Lat, color=factor(Year)), position = position_jitter(width = 0.01, height = 0.01), alpha = 0.7, size = .02) +
#   scale_color_manual(values = c("2015" = "red", "1990" = "black")) +
#   theme_void() + theme(legend.position="none")
# 
# #Save plot
# pdf("bt_isolate_map_inset.pdf", height = 5, width = 9)
# grid.newpage()
# v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
# v2 <- viewport(width = 0.1, height = 0.1, x = 0.5, y = 0.9) #plot area for the inset map
# print(map, vp=v1) 
# print(ramsey_map, vp=v2)
# dev.off()

#Previous analysis with FAO data
#Oat data was downloaded from FAOSTAT on June 25th, 2018. All available data (area harvested, yield, and production) was downloaded for 1990 - 2015.
# oat_data <- read.csv("oat_raw_data/FAOSTAT_data_6-25-2018.csv", header = TRUE)
# 
# #Set theme to be classic
# theme_set(theme_classic() + 
#             theme(
#               legend.text=element_text(size=11),
#               axis.text=element_text(size=9),
#               axis.title = element_text(size=12),
#               axis.line = element_line(size=0.5)
#               #strip.text.y = element_blank()
#             )
# )
# 
# oat_plot <- ggplot(oat_data %>% filter(Element != "Yield"), aes(x = Year, y = Value)) +
#   geom_line() + 
#   scale_y_continuous(label=comma) +
#   facet_wrap(~Element, strip.position = "left", scales = "free_y", nrow = 2, labeller = as_labeller(c("Area harvested" = "Area Harvested (ha)", "Production" = "Production (tonnes)"))) +
#   ylab(NULL) +
#   theme(strip.background = element_blank(), strip.placement = "outside", strip.text.y = element_text(size = 12), panel.spacing = unit(1.5, "lines"))
# 
# #Save plot
# pdf("oat_production.pdf", height = 4, width = 6)
# oat_plot
# dev.off()
 