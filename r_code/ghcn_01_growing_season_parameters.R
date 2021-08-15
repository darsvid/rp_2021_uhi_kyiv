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
    #here we select the longest gaps only for March-October
    #since long growing season lays within these time boundaries
    filter(month(start) %in% c(3:10) & gap_days >= 5) 

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
#no need to introduce |is.na(TAVG_ard) because we filled all the NA gaps with linear
#interpolation
long_dates_df <- long_yes_df[,-2:-3]
long_dates_df <- rbind(c(NA, ifelse(long_yes_df$LONG_GS[1] == 1, yes = 0, no = NA)),
                       long_dates_df,
                       c(NA, ifelse(long_yes_df$LONG_GS[length(long_yes_df$LONG_GS)] == 1, 0, NA)))

long_periods_df <- data.frame(
    start = long_dates_df$DATE[which(diff(long_dates_df$LONG_GS, lag = 1) == 1) + 1],
    end =   long_dates_df$DATE[which(diff(long_dates_df$LONG_GS, lag = 1) == -1)])
    
long_gs_df <- long_periods_df %>%
    mutate(length_days = as.numeric(end - start) + 1) %>%
    group_by(year = year(start)) %>%
    filter(length_days == max(length_days)) %>%
    mutate(start_doy = as.numeric(yday(start)), end_doy = as.numeric(yday(end))) %>%
    select(year, start, end, start_doy, end_doy, length_days)

saveRDS(long_gs_df, paste0(output_dir, "long_gs_df.rds"))

# active growing season (>=10)
active_yes_df <- ard_df %>%
    mutate(ACTIVE_GS = ifelse(TAVG_ard < 10, 0, 1))
#no need to introduce |is.na(TAVG_ard) because we filled all the NA gaps with linear
#interpolation
active_dates_df <- active_yes_df[,-2:-3]
active_dates_df <- rbind(c(NA, ifelse(active_yes_df$ACTIVE_GS[1] == 1, yes = 0, no = NA)),
                       active_dates_df,
                       c(NA, ifelse(active_yes_df$ACTIVE_GS[length(active_yes_df$ACTIVE_GS)] == 1, 0, NA)))

active_periods_df <- data.frame(
    start = active_dates_df$DATE[which(diff(active_yes_df$ACTIVE_GS, lag = 1) == 1) + 1],
    end =   active_dates_df$DATE[which(diff(active_yes_df$ACTIVE_GS, lag = 1) == -1)])

active_gs_df <- active_periods_df %>%
    mutate(length_days = as.numeric(end - start) + 1) %>%
    group_by(year = year(start)) %>%
    filter(length_days == max(length_days)) %>%
    mutate(start_doy = as.numeric(yday(start)), end_doy = as.numeric(yday(end))) %>%
    select(year, start, end, start_doy, end_doy, length_days)

saveRDS(active_gs_df, paste0(output_dir, "active_gs_df.rds"))

# summer growing season (>=15)
summer_yes_df <- ard_df %>%
    mutate(SUMMER_GS = ifelse(TAVG_ard < 15, 0, 1))
#no need to introduce |is.na(TAVG_ard) because we filled all the NA gaps with linear
#interpolation
summer_dates_df <- summer_yes_df[,-2:-3]
summer_dates_df <- rbind(c(NA, ifelse(summer_yes_df$SUMMER_GS[1] == 1, yes = 0, no = NA)),
                         summer_dates_df,
                         c(NA, ifelse(summer_yes_df$SUMMER_GS[length(summer_yes_df$SUMMER_GS)] == 1, 0, NA)))

summer_periods_df <- data.frame(
    start = summer_dates_df$DATE[which(diff(summer_yes_df$SUMMER_GS, lag = 1) == 1) + 1],
    end =   summer_dates_df$DATE[which(diff(summer_yes_df$SUMMER_GS, lag = 1) == -1)])

summer_gs_df <- summer_periods_df %>%
    mutate(length_days = as.numeric(end - start) + 1) %>%
    group_by(year = year(start)) %>%
    filter(length_days == max(length_days)) %>%
    mutate(start_doy = as.numeric(yday(start)), end_doy = as.numeric(yday(end))) %>%
    select(year, start, end, start_doy, end_doy, length_days)

saveRDS(summer_gs_df, paste0(output_dir, "summer_gs_df.rds"))