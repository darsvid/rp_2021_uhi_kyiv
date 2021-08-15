library(tidyverse)
library(mblm)

input_dir  <- "F:/rp_2021_uhi_kyiv/data_output/GHCN_analysis_results/"
setwd(input_dir)

long_gs_df   <- readRDS("long_gs_df.rds")
active_gs_df <- readRDS("active_gs_df.rds")
summer_gs_df <- readRDS("summer_gs_df.rds") %>% group_by(year) %>%
    filter(row_number()==1) #to remove ties 1956, 1979, 1983, 1984


col_names  <- c("year", "long_start", "long_end", "long_length",
                "active_start", "active_end", "active_length",
                "summer_start", "summer_end", "summer_length")

all_gs_df <- cbind(long_gs_df[,c(-2:-3)], active_gs_df[,c(-1:-3)], summer_gs_df[,c(-1:-3)])
colnames(all_gs_df) <- col_names

sen <- function(..., weights = NULL) {
    mblm::mblm(..., repeated = F)
}

sieg <- function(..., weights = NULL) {
    mblm::mblm(...)
}

#----- start of growing season plot -------
start_long_df <- all_gs_df %>% select(year, long_start, active_start, summer_start) %>%
    pivot_longer(cols = -year, names_to = "growing_season", values_to = "day") %>%
    mutate(parameter = factor(growing_season,
                              levels = c("long_start", "active_start", "summer_start")),
           period = "1881-2020") %>%
    select(period, year, parameter, day)

period_02 <- start_long_df %>%
    filter(year %in% 1961:2020) %>%
    mutate(period = "1961-2020")

period_03 <- start_long_df %>%
    filter(year %in% 1961:1990) %>%
    mutate(period = "1961-1990")

period_04 <- start_long_df %>%
    filter(year %in% 1991:2020) %>%
    mutate(period = "1991-2020")

start_by_period_df <- rbind(start_long_df, period_02, period_03, period_04) %>%
    mutate(periods = factor(period,
                            levels = c("1881-2020", "1961-2020", "1961-1990", "1991-2020")),
           combine = paste0(periods,"-", parameter),
           period_param = factor(combine,
                                 levels = c("1881-2020-long_start", "1881-2020-active_start", "1881-2020-summer_start",
                                            "1961-2020-long_start", "1961-2020-active_start", "1961-2020-summer_start",
                                            "1961-1990-long_start", "1961-1990-active_start", "1961-1990-summer_start",
                                            "1991-2020-long_start", "1991-2020-active_start", "1991-2020-summer_start"))) %>%
    select(period_param, periods, parameter, year, day)

param_names <- c(
    "1881-2020-long_start" = "T \u2265 5°C",
    "1881-2020-active_start" = "T \u2265 10°C",
    "1881-2020-summer_start" = "T \u2265 15°C",
    "1961-2020-long_start" = "T \u2265 5°C",
    "1961-2020-active_start" = "T \u2265 10°C",
    "1961-2020-summer_start" = "T \u2265 15°C",
    "1961-1990-long_start" = "T \u2265 5°C",
    "1961-1990-active_start" = "T \u2265 10°C",
    "1961-1990-summer_start" = "T \u2265 15°C",
    "1991-2020-long_start" = "T \u2265 5°C",
    "1991-2020-active_start" = "T \u2265 10°C",
    "1991-2020-summer_start" = "T \u2265 15°C"
)


png("F:/rp_2021_uhi_kyiv/fig_output/grow_start_param.png", width=16, height=16,
    units = "cm", res = 500)

plot_start <- ggplot(start_by_period_df, aes(year, day)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_smooth(method = "lm", se = F, aes(colour = "blue")) +
    geom_smooth(method = sen, se = F, aes(colour = "orange"), linetype = "longdash") +
    geom_smooth(method = sieg, se = F, aes(colour = "red"), linetype = "longdash") +
    facet_wrap(~period_param, ncol=3, scales = "free", labeller = as_labeller(param_names)) +
    theme_bw() +
    labs(y = "day of year") +
    theme(axis.text=element_text(size=8, colour = "black"), axis.title.x = element_blank(),
          #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 10, colour = "black"),
          legend.position = "bottom") +
    scale_color_identity(guide = "legend",
                         name = "Estimator",
                         breaks = c("blue", "orange", "red"),
                         labels = c("Ordinary least squares", "Theil-Sen single median",
                                    "Siegel repeated medians"))
plot_start

dev.off()


#-----end of growing season plot------
end_long_df <- all_gs_df %>% select(year, long_end, active_end, summer_end) %>%
    pivot_longer(cols = -year, names_to = "growing_season", values_to = "day") %>%
    mutate(parameter = factor(growing_season,
                              levels = c("long_end", "active_end", "summer_end")),
           period = "1881-2020") %>%
    select(period, year, parameter, day)

period_02 <- end_long_df %>%
    filter(year %in% 1961:2020) %>%
    mutate(period = "1961-2020")

period_03 <- end_long_df %>%
    filter(year %in% 1961:1990) %>%
    mutate(period = "1961-1990")

period_04 <- end_long_df %>%
    filter(year %in% 1991:2020) %>%
    mutate(period = "1991-2020")

end_by_period_df <- rbind(end_long_df, period_02, period_03, period_04) %>%
    mutate(periods = factor(period,
                            levels = c("1881-2020", "1961-2020", "1961-1990", "1991-2020")),
           combine = paste0(periods,"-", parameter),
           period_param = factor(combine,
                                 levels = c("1881-2020-long_end", "1881-2020-active_end", "1881-2020-summer_end",
                                            "1961-2020-long_end", "1961-2020-active_end", "1961-2020-summer_end",
                                            "1961-1990-long_end", "1961-1990-active_end", "1961-1990-summer_end",
                                            "1991-2020-long_end", "1991-2020-active_end", "1991-2020-summer_end"))) %>%
    select(period_param, periods, parameter, year, day)

param_names <- c(
    "1881-2020-long_end" = "T \u2265 5°C",
    "1881-2020-active_end" = "T \u2265 10°C",
    "1881-2020-summer_end" = "T \u2265 15°C",
    "1961-2020-long_end" = "T \u2265 5°C",
    "1961-2020-active_end" = "T \u2265 10°C",
    "1961-2020-summer_end" = "T \u2265 15°C",
    "1961-1990-long_end" = "T \u2265 5°C",
    "1961-1990-active_end" = "T \u2265 10°C",
    "1961-1990-summer_end" = "T \u2265 15°C",
    "1991-2020-long_end" = "T \u2265 5°C",
    "1991-2020-active_end" = "T \u2265 10°C",
    "1991-2020-summer_end" = "T \u2265 15°C"
)

png("F:/rp_2021_uhi_kyiv/fig_output/grow_end_param.png", width=16, height=16,
    units = "cm", res = 500)

plot_end <- ggplot(end_by_period_df, aes(year, day)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_smooth(method = "lm", se = F, aes(colour = "blue")) +
    geom_smooth(method = sen, se = F, aes(colour = "orange"), linetype = "longdash") +
    geom_smooth(method = sieg, se = F, aes(colour = "red"), linetype = "longdash") +
    facet_wrap(~period_param, ncol=3, scales = "free", labeller = as_labeller(param_names)) +
    theme_bw() +
    labs(y = "day of year") +
    theme(axis.text=element_text(size=8, colour = "black"), axis.title.x = element_blank(),
          #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 10, colour = "black"),
          legend.position = "bottom") +
    scale_color_identity(guide = "legend",
                         name = "Estimator",
                         breaks = c("blue", "orange", "red"),
                         labels = c("Ordinary least squares", "Theil-Sen single median",
                                    "Siegel repeated medians")) #face = "bold"
plot_end

dev.off()


#-----length of growing season plot------
length_long_df <- all_gs_df %>% select(year, long_length, active_length, summer_length) %>%
    pivot_longer(cols = -year, names_to = "growing_season", values_to = "day") %>%
    mutate(parameter = factor(growing_season,
                              levels = c("long_length", "active_length", "summer_length")),
           period = "1881-2020") %>%
    select(period, year, parameter, day)

period_02 <- length_long_df %>%
    filter(year %in% 1961:2020) %>%
    mutate(period = "1961-2020")

period_03 <- length_long_df %>%
    filter(year %in% 1961:1990) %>%
    mutate(period = "1961-1990")

period_04 <- length_long_df %>%
    filter(year %in% 1991:2020) %>%
    mutate(period = "1991-2020")

length_by_period_df <- rbind(length_long_df, period_02, period_03, period_04) %>%
    mutate(periods = factor(period,
                            levels = c("1881-2020", "1961-2020", "1961-1990", "1991-2020")),
           combine = paste0(periods,"-", parameter),
           period_param = factor(combine,
                                 levels = c("1881-2020-long_length", "1881-2020-active_length", "1881-2020-summer_length",
                                            "1961-2020-long_length", "1961-2020-active_length", "1961-2020-summer_length",
                                            "1961-1990-long_length", "1961-1990-active_length", "1961-1990-summer_length",
                                            "1991-2020-long_length", "1991-2020-active_length", "1991-2020-summer_length"))) %>%
    select(period_param, periods, parameter, year, day)

param_names <- c(
    "1881-2020-long_length" = "T \u2265 5°C",
    "1881-2020-active_length" = "T \u2265 10°C",
    "1881-2020-summer_length" = "T \u2265 15°C",
    "1961-2020-long_length" = "T \u2265 5°C",
    "1961-2020-active_length" = "T \u2265 10°C",
    "1961-2020-summer_length" = "T \u2265 15°C",
    "1961-1990-long_length" = "T \u2265 5°C",
    "1961-1990-active_length" = "T \u2265 10°C",
    "1961-1990-summer_length" = "T \u2265 15°C",
    "1991-2020-long_length" = "T \u2265 5°C",
    "1991-2020-active_length" = "T \u2265 10°C",
    "1991-2020-summer_length" = "T \u2265 15°C"
)

png("F:/rp_2021_uhi_kyiv/fig_output/grow_length_param.png", width=16, height=16,
    units = "cm", res = 500)

plot_end <- ggplot(length_by_period_df, aes(year, day)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_smooth(method = "lm", se = F, aes(colour = "blue")) +
    geom_smooth(method = sen, se = F, aes(colour = "orange"), linetype = "longdash") +
    geom_smooth(method = sieg, se = F, aes(colour = "red"), linetype = "longdash") +
    facet_wrap(~period_param, ncol=3, scales = "free", labeller = as_labeller(param_names)) +
    theme_bw() +
    labs(y = "days") +
    theme(axis.text=element_text(size=8, colour = "black"), axis.title.x = element_blank(),
          #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 10, colour = "black"),
          legend.position = "bottom") +
    scale_color_identity(guide = "legend",
                         name = "Estimator",
                         breaks = c("blue", "orange", "red"),
                         labels = c("Ordinary least squares", "Theil-Sen single median",
                                    "Siegel repeated medians")) #face = "bold"
plot_end

dev.off()

# plot <- ggplot(all_gs_df, aes(year)) +
#     geom_ribbon(aes(ymin =long_start, ymax = long_end), fill = "grey50") +
#     geom_ribbon(aes(ymin =active_start, ymax = active_end), fill = "green") +
#     geom_ribbon(aes(ymin =summer_start, ymax = summer_end), fill = "red") +
#     theme_bw()
# 
# 
# all_gs_long_df <- all_gs_df%>%
#     pivot_longer(cols = -year, names_to = "growing_season", values_to = "day") #%>%
#mutate(param = factor(growing_season, levels = c("start_doy", "end_doy", "length")) )