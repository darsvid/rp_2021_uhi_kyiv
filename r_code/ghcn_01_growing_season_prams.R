library(tidyverse)
library(lubridate)
library(modern)
library(zoo)

input_dir  <- "F:/rp_2021_uhi_kyiv/data_raw/GHCN_climate_record"
output_dir <- "F:/rp_2021_uhi_kyiv/data_output/GHCN_analysis_results/"


setwd(input_dir)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
}

raw_df       <- read.csv("ghcd_kyiv_1881-2020.csv", colClasses=c("DATE"="Date"),
                         stringsAsFactors = F) %>%
                select(DATE, TAVG) %>%
                mutate(TAVG_noout = iglewicz_hoaglin(TAVG))
all_dates_df <- as.data.frame(seq.Date(min(raw_df$DATE), max(raw_df$DATE),
                                      by = "day"))
colnames(all_dates_df) <- "DATE"


grow_df <- full_join(all_dates_df, raw_df) %>%
    #filter(month(DATE) %in% c(3:10)) %>%
    select(DATE, TAVG_noout) %>%
    mutate(GAP = ifelse(is.na(TAVG_noout), 1, 0))

# if the number of consecutive gap days in a certain month >= 5 then
# exclude this year from the analysis
gap_df <- grow_df[,-2]
gap_df <- rbind(c(NA, ifelse(gap_df$GAP[1] == 1, yes = 0, no = NA)),
                gap_df,
                c(NA, ifelse(gap_df$GAP[gap_df$GAP] == 1, 0, NA)))
    
#
gap_periods_df <- data.frame(
    start = gap_df$DATE[which(diff(gap_df$GAP, lag = 1) == 1) + 1],
    end = gap_df$DATE[which(diff(gap_df$GAP, lag = 1) == -1)]) %>%
    mutate(gap_days = as.numeric(end - start) + 1) %>%
    filter(gap_days >= 5)

excl_years <- unique(year(gap_periods_df$start)) 

ard_df <- grow_df[,-3] %>%
    filter(!year(DATE) %in% excl_years) %>%
    mutate(TAVG_fill = na.approx(TAVG_noout)) %>% # fill small gaps with linear interpolation
    rename(TAVG_ard = TAVG_fill)

# important - you may still analyze longest gaps for months 3-10 when it comes to exclude years
# or caculate them by months

#define the parameters of growing season
# long (>=5), active (>=10) and summer (>=15)

# long growing season (>=5)
long_yes_df <- ard_df %>%
    mutate(LONG_GS = ifelse(TAVG_ard < 5, 0, 1))
#no need to introduce |is.na(TAVG) because we filled all the NA gaps with linear
#interpolation

long_dates_df <- rbind(c(NA, ifelse(long_yes_df$LONG_GS[1] == 1, yes = 0, no = NA)),
                       long_yes_df[,-2],
                       c(NA, ifelse(long_yes_df$LONG_GS[length(long_yes_df$LONG_GS)] == 1, 0, NA)))
long_periods_df <- data.frame(
    start = long_dates_df$DATE[which(diff(long_dates_df$LONG_GS, lag = 1) == 1) + 1],
    end = long_dates_df$DATE[which(diff(long_dates_df$LONG_GS, lag = 1) == -1)])

long_gs_df <- long_periods_df %>%
    mutate(length_days = as.numeric(end - start) + 1) %>%
    group_by(year = year(start)) %>%
    filter(length_days == max(length_days)) %>%
    mutate(start_doy = as.numeric(yday(start)), end_doy = as.numeric(yday(end))) %>%
    select(year, start, end, start_doy, end_doy, length_days)
    #rename()



# compl_month_df <- ghcn_modern_df %>%
#     group_by(year = year(DATE), month = month(DATE), days = days_in_month(DATE)) %>%
#     summarise(obs_n = n(), na_n = sum(is.na(TAVG))) %>%
#     mutate(empty_days = days - (obs_n - na_n))
# exclude_years <- unique(compl_month_df$year[compl_month_df$empty_days<5])
# 
# ard_modern_df <- ghcn_modern_df %>%
#     filter(!year(DATE) %in% exclude_years)
    