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

linear_df <- data.frame(period = character(), parameter = character(),
                        slope = numeric(), r_sq = numeric(), r_sq_adj = numeric(),
                        p_value = numeric(), ts_slope = numeric(), ts_p = numeric(),
                        sgl_slope = numeric(), sgl_p = numeric())
periods <- c("1881:2020", "1961:2020", "1961:1990", "1991:2020")

for (i in 2:ncol(all_gs_df)){
    for (j in 1:length(periods)){
        # define time period
        p <- strsplit(periods[j], ":")
        t <- as.integer(p[[1]][1]):as.integer(p[[1]][2])
        #define variables
        temp_df <- all_gs_df[,c(1,i)] %>%
            filter(year %in% t)
        x <- temp_df$year
        y <- temp_df[,2][[1]]
        # calculate different linear models
        # linear model and its parameters
        linear <- lm(y ~ x)
        slope <- linear$coefficients[[2]]
        r_sq  <- summary(linear)$r.squared
        r_sq_adj <- summary(linear)$adj.r.squared
        p_value <-  summary(linear)$coefficients[,4][2]
        # theil-sen model single median
        ts <- mblm(y ~ x, repeated = F)
        ts_slope <- ts$coefficients[[2]]
        ts_p     <- summary(ts)$coefficients[,4][2]
        # siegel model repeated medians
        sgl <- mblm(y ~ x, repeated = T)
        sgl_slope <- sgl$coefficients[[2]]
        sgl_p     <- summary(sgl)$coefficients[,4][2]
        
        params <- c(periods[j], colnames(all_gs_df[,i]), slope, r_sq, r_sq_adj, p_value,
                    ts_slope, ts_p, sgl_slope, sgl_p)
        linear_df <- rbind(linear_df, params)
    }
       
}

colnames(linear_df) <- c("period", "parameter", "slope", "r_sq", "r_sq_adj", "p_value",
                         "ts_slope", "ts_p", "sgl_slope", "sgl_p")
saveRDS(linear_df,   "linear_trends_df.rds")
write.csv(linear_df, "linear_trends_df.csv")


# predict active growing season changes
# using on Theil-Sen single median slope model

fit_df    <- data.frame(model = character(), parameter = character(),
                        fit_year = numeric(), fit = numeric(), lwr = numeric(),
                        upr = numeric())
periods   <- c("1881:2020", "1961:2020", "1991:2020")
fit_years <- data.frame(x = c(2021, 2030, 2040, 2050))


for (i in 2:ncol(all_gs_df)){
    for (j in 1:length(periods)){
        # define time period
        p <- strsplit(periods[j], ":")
        t <- as.integer(p[[1]][1]):as.integer(p[[1]][2])
        #define variables
        temp_df <- all_gs_df[,c(1,i)] %>%
            filter(year %in% t)
        x <- temp_df$year
        y <- temp_df[,2][[1]]
        # Theil-Sen model single median
        ts <- mblm(y ~ x, repeated = F)
        pred_df <- predict(ts, newdata = fit_years, interval = "confidence")
        
        param <- colnames(all_gs_df[,i])
        
        temp_df <- cbind(rep(periods[j], nrow(fit_years)), rep(param, nrow(fit_years)),
                         fit_years, pred_df)
        
        fit_df <- rbind(fit_df, temp_df)
    }
    
}

colnames(fit_df) <- c("period", "parameter", "fit_year", "fit", "lwr", "upr")
saveRDS(fit_df,   "ts_model_prediction_df.rds")
write.csv(fit_df, "ts_model_prediction_df.csv")