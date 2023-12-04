"""
Calibration.

Script to calibrate things. Contains reference parameterization.

"""
import os
# import csv
import configparser
# import pandas as pd
# import geopandas
# from collections import OrderedDict
import numpy as np
import pandas as pd 

from misc import shell_volume

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')

### Construct the elemental parameter list
load_params_element = {
	###
	'k': 25000, # people-mb/sat-s
	'lambda': 2, # ground station latency coefficient. unitless
	'mu': 30, # unavoidable minimum latency coefficient. ms
	'atmo_damages': 0,
	'a': 9,
	'a_L': 6, # $/y/ms # 
	'a_S': 150, #153.238 # $/Mb/s
	'a_U': 2, # preference parameter for uptime
	# a_L and a_S are from Liu et al 2018
	'L_bar': 275, #550/3, # maximum acceptable latency
	'S_bar': 55, #100.25, # testing
	# calibrated from existing plans terrestrially
	'theta_ubar': 1/2, # unitless demand location parameter
	'theta_star': (1 + 2*1/2)/3, # cutoff consumer utility level, (1 + 2*theta_ubar)/3
	'tau': 2, #1/5, # seconds/maneuver.
	'eta': 1, # initialization so that the congestion function works
	'N': 1e7, # people . 100 million people in the mature market.
	'c': 400000, #1.5e5, # $/sat. constant term in unit cost of launching sat
	'd': 1000, #300, # $/sat/km . linear term in unit cost of launching sat
	'e': 2, #1, # $/sat/km^2 . quadratic term in unit cost of launching sat
	'K': 1, # scaling factor for costs
	# d, e calibrated to make 500 km the unit-cost-minimizing altitude
	# c calibrated below given d, e to make 5-year cost of a satellite equal $500,000. Annualized to one year with no discounting, we get $100,000.
	'radius': 1, #1, # km. Safety margin for collision avoidance maneuvers. Calibrated.
	'mass': 223, # kg. mass of satellites. From Lifson et al 2022. Not currently being used anywhere, but could be used to calculate cost/kg.
	'rho': 17.5, # km. shell thickness
	'objects_sh': 0 #NULL #objects_distns, # initialize the objects distribution
}

# Hugh Lewis analysis: 
# https://twitter.com/ProfHughLewis/status/1478316593675395072. ~20 maneuvers per day in year ending 2021
# https://twitter.com/ProfHughLewis/status/1610307514792382471. ~75 maneuvers per day in year ending 2022.

# ### Check (operational) Starlink altitude distribution
# read_json("../data/leo-satcat-26-12-2022.json", simplifyVector=TRUE) %>%
# 	filter(CURRENT=="Y", str_detect(OBJECT_NAME,"STARLINK"), PERIGEE > 500) %>%
# 	select(APOGEE, PERIGEE) %>% 
# 	mutate(mean_alt = (as.numeric(APOGEE) + as.numeric(PERIGEE))/2) %>% 
# 	select(mean_alt) %>%
# 	summary()

### Read in jspoc object count data
# object_counts_2022 <- read_json("../data/leo-satcat-26-12-2022.json", simplifyVector=TRUE) %>% 
# 	filter(CURRENT=="Y" & APOGEE < 550+35/2 & PERIGEE > 550-35/2)
# total_count <- nrow(object_counts_2022)
# starlink_count <- object_counts_2022 %>% filter(str_detect(OBJECT_NAME,"STARLINK")) %>% nrow()
# total_count
# starlink_count

def get_object_counts():
    """
    Hugh Lewis analysis: 
    https://twitter.com/ProfHughLewis/status/1478316593675395072. ~20 maneuvers per day in year ending 2021
    https://twitter.com/ProfHughLewis/status/1610307514792382471. ~75 maneuvers per day in year ending 2022.

    original_paper: 

    total_count = 5616
    starlink_count = 3015

    this function:

    total_count = 3149 
    starlink_count = 3015

    """
    object_counts_2022 = pd.read_json(os.path.join(DATA_RAW, 'leo-satcat-26-12-2022.json'))
    if not os.path.exists(os.path.join(DATA_PROCESSED, 'leo-satcat-26-12-2022.csv')):
        object_counts_2022.to_csv(os.path.join(DATA_PROCESSED, 'leo-satcat-26-12-2022.csv'), index=False)

    object_counts_2022 = object_counts_2022.to_dict('records')

    upper_bound = 550+35/2
    lower_bound = 550-35/2
    
    total = []
    starlink_count = 0

    for item in object_counts_2022:

        if item['CURRENT'] == 'Y' and float(item['PERIGEE']) > lower_bound and float(item['APOGEE']) < upper_bound:

            total.append({
                'INTLDES': item['INTLDES'],
                'NORAD_CAT_ID': item['NORAD_CAT_ID'],
                'OBJECT_NAME': item['OBJECT_NAME'],
                'APOGEE': item['APOGEE'],
                'PERIGEE': item['PERIGEE'],
            })

            if 'STARLINK' in item['OBJECT_NAME']:
                starlink_count += 1 

    total = pd.DataFrame(total)
    total.to_csv(os.path.join(DATA_PROCESSED, 'leo-satcat-26-12-2022_selected.csv'), index=False)

    total_count = len(total)

    return total_count, starlink_count


# params_element$radius <- sqrt((75 * shell_volume(550)) / ( pi *2 * 86400 * total_count * starlink_count * velocity(550)))
# params_element$radius

# ### Produce simplified CSV of object counts: 3 groups, Starlink / Other / large objects
# object_sheet <- read_json("../data/leo-satcat-26-12-2022.json", simplifyVector=TRUE) |> #nolint
# 	filter(CURRENT=="Y") |>
# 	mutate(mean_altitude = (as.numeric(APOGEE)+as.numeric(PERIGEE))/2) |>
# 	select(OBJECT_TYPE, OBJECT_NAME, RCS_SIZE, mean_altitude) |>
# 	mutate(
# 		OBJECT_TYPE = ifelse(str_detect(OBJECT_NAME,"STARLINK"), "STARLINK", OBJECT_TYPE),
# 		OBJECT_TYPE = ifelse(str_detect(OBJECT_NAME,"ONEWEB"), "ONEWEB", OBJECT_TYPE)
# 	)

# summary(object_sheet$OBJECT_TYPE |> as.factor())

# # Create a new dataset with object count for each group based on altitude bins and object type. _all is what's used for simulations, _nosmall is for comparison.
# binned_data_nosmall <- object_sheet %>%
# 	filter(RCS_SIZE == "MEDIUM" | RCS_SIZE == "LARGE") %>%
# 	mutate(altitude_bin = cut(mean_altitude, breaks = seq(200, max(mean_altitude), by = 35), include.lowest = TRUE)) %>%
# 	group_by(altitude_bin, OBJECT_TYPE) %>%
# 	summarise(object_count = n(), .groups = 'drop') %>%
# 	spread(key = OBJECT_TYPE, value = object_count, fill = 0) %>%
# 	mutate(lower_bound = as.numeric(str_extract(as.character(altitude_bin), "[-+]?\\d*\\.?\\d+([eE][-+]?\\d+)?(?=,)")),
#          upper_bound = as.numeric(str_extract(as.character(altitude_bin), "(?<=,)[-+]?\\d*\\.?\\d+([eE][-+]?\\d+)?")),
#          bin_center = (lower_bound + upper_bound) / 2) %>%
# 	select(-altitude_bin, Altitude_Bin = altitude_bin, Bin_Center = bin_center, everything()) |>
# 	filter(!is.na(Bin_Center))

# binned_data_all <- object_sheet %>%
# 	mutate(altitude_bin = cut(mean_altitude, breaks = seq(200, max(mean_altitude), by = (params_element$rho)), include.lowest = TRUE)) %>%
# 	group_by(altitude_bin, OBJECT_TYPE) %>%
# 	summarise(object_count = n(), .groups = 'drop') %>%
# 	spread(key = OBJECT_TYPE, value = object_count, fill = 0) %>%
# 	mutate(lower_bound = as.numeric(str_extract(as.character(altitude_bin), "[-+]?\\d*\\.?\\d+([eE][-+]?\\d+)?(?=,)")),
#          upper_bound = as.numeric(str_extract(as.character(altitude_bin), "(?<=,)[-+]?\\d*\\.?\\d+([eE][-+]?\\d+)?")),
#          bin_center = (lower_bound + upper_bound) / 2) %>%
# 	select(-altitude_bin, Altitude_Bin = altitude_bin, Bin_Center = bin_center, everything()) |>
# 	filter(!is.na(Bin_Center))

# write_csv(binned_data_nosmall, "../data/binned_counts_nosmall.csv")
# write_csv(binned_data_all, "../data/binned_counts_all.csv")

# ### Construct the objects_sh object from binned_data_all. The objects_sh object has two columns: h and objects. h is Bin_Center, and objects is the row-wise sum of DEBRIS, PAYLOAD, `ROCKET BODY`, and UNKNOWN. Column names are h and objects.
# objects_sh <- binned_data_all %>%
#   mutate(objects = DEBRIS + PAYLOAD + `ROCKET BODY` + UNKNOWN) %>%
#   select(Bin_Center, objects) %>%
#   rename(h = Bin_Center)

# params_element$objects_sh <- objects_sh # Update the params_element object with the objects_sh object

# ### Make a stacked bar chart of object counts in binned_data_all using ggplot. The stacking should be done using the DEBRIS, PAYLOAD, `ROCKET BODY`, STARLINK, and UNKNOWN variables.
# all_histogram <- binned_data_all %>%
#   gather(key = "object_type", value = "count", DEBRIS, PAYLOAD, `ROCKET BODY`, STARLINK, ONEWEB, UNKNOWN) %>%
#   ggplot(aes(x = Bin_Center, y = count, fill = object_type)) +
#   geom_bar(stat = "identity", width = 35, position = "stack") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   xlab("Altitude (km)") +
#   ylab("Object Count") +
#   ggtitle("All tracked objects") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# nosmall_histogram <- binned_data_nosmall %>%
#   gather(key = "object_type", value = "count", DEBRIS, PAYLOAD, `ROCKET BODY`, STARLINK, ONEWEB, UNKNOWN) %>%
#   ggplot(aes(x = Bin_Center, y = count, fill = object_type)) +
#   geom_bar(stat = "identity", width = 35, position = "stack") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   xlab("Altitude (km)") +
#   ylab("Object Count") +
#   ggtitle("Only MEDIUM or LARGE objects (RCS >= 0.1 m^2)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave(
# 	"../figures/object_count_histograms.png", 
# 	all_histogram / nosmall_histogram, 
# 	width = 10, 
# 	height = 10,
# 	bg="white")

if __name__ == '__main__':

    # params_element = load_params_element()

    total_count, starlink_count = get_object_counts()
    print(total_count, starlink_count)

    
    # params_element['radius'] = np.sqrt(
    #     (75 * shell_volume(550)) / 
    #     ( np.pi *2 * 86400 * total_count * starlink_count * velocity(550))
    #     )

