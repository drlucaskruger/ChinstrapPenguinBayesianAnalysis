

#  Load libraries -------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(brms)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(stringr)
library(purrr)
library(cluster)
library(factoextra)
library(ggplotify)

# -----------Read and prepare data ------------------------------------------------
all <- read.csv("C:/MAPPPD/AllCounts_V_4_4.csv") %>% 
  filter(season_starting >= 1970)

GENERATION_LENGTH <- 9.4   # chinstrap penguin generation length (years)

# Initial data preparation (including zeros, accuracy <=5)
model_data <- all %>%
  filter(common_name == "chinstrap penguin",
         count_type == "nests",
         accuracy <= 5) %>%        # keep all accuracy levels, zeros included
                                  # 
  group_by(site_id) %>%
  mutate(
    mean_site_size = mean(penguin_count, na.rm = TRUE),
    centered_year = season_starting - mean(season_starting),
    log_count = log(penguin_count + 1),   # log(0+1)=0
    site_id = as.factor(site_id),
    accuracy_factor = as.factor(accuracy)
  ) %>%
  group_by(site_id) %>%
  filter(n_distinct(season_starting) >= 2) %>%   # at least 2 years per site
  
  ungroup()%>%
  filter(cammlr_region %in% c("48.1","48.2"))

# --------- Identify and exclude problematic sites/observations -----------------


# Based on manual inspection of accuracy 4/5 data
problematic_sites <- c("GOWC", "KELL", "OBRI", "TART")
exclusion_log <- data.frame(
  site_id = character(),
  year = integer(),
  reason = character(),
  action = character(),
  stringsAsFactors = FALSE
)

# Remove entire problematic sites
model_data_clean <- model_data %>%
  filter(!site_id %in% problematic_sites)

# Remove single extreme observation: Harmony Point 1971
model_data_clean <- model_data_clean %>%
  filter(!(site_id == "HARM" & season_starting == 1971))





# 
message("Original observations: ", nrow(model_data))
message("Cleaned observations: ", nrow(model_data_clean))

# Hierarchical clustering to define size categories -------------------
colony_cluster_data <- model_data_clean %>%
  group_by(site_id) %>%
  summarise(mean_size = mean(penguin_count, na.rm = TRUE)) %>%
  mutate(log_size = log10(mean_size + 1),
         scaled_size = scale(log_size))

# Hierarchical clustering with complete linkage



# Compute the gap statistic for up to 10 possible clusters
set.seed(123)
gap_stat <- clusGap((colony_cluster_data$scaled_size),
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 500)

# Plot the gap statistic
plot(gap_stat)


#save("gap_analysis.pdf",width = 6, height = 4, dpi = 600) # doesn't work 


hc <- hclust(dist(colony_cluster_data$log_size), method = "complete")


# Choose k=4 based on gap statistic 


plot(hc)

##ggsave("cluster.pdf", width = 6, height = 4, dpi = 600)

k_clusters <- 4

colony_cluster_data$size_cluster <- cutree(hc, k = k_clusters)

# Map cluster numbers to intuitive names (ordered by increasing size)
cluster_means <- colony_cluster_data %>%
  group_by(size_cluster) %>%
  summarise(mean_size = mean(mean_size)) %>%
  arrange(mean_size) %>%
  mutate(size_category = c("tiny", "small", "typical", "large"))
                          
print(cluster_means)

# Join cluster names back to main data
colony_cluster_data <- colony_cluster_data %>%
  left_join(cluster_means[, c("size_cluster", "size_category")], by = "size_cluster")

ggplot(colony_cluster_data,aes(reorder(size_category,-mean_size),mean_size))+
  geom_boxplot()+
  scale_y_log10(labels=scales::comma_format())+
  theme_bw()+
  xlab("Colony size category")+
  ylab("Mean number of active nests")
  
##ggsave("colony_size_categories.pdf", width = 6, height = 4, dpi = 600)

colony_size_ranges<-colony_cluster_data%>%
  group_by(size_category)%>%
  summarise(min=min(mean_size),
            max=max(mean_size))

print(colony_size_ranges)

table(colony_cluster_data$size_cluster)


model_data_clean <- model_data_clean %>%
  left_join(colony_cluster_data[, c("site_id", "size_category")], by = "site_id")


table(model_data_clean$size_category)


# Ensure factor order
model_data_clean$size_category <- factor(model_data_clean$size_category,
                                         levels = c("tiny", "small", "typical", "large"))

# 5. Final dataset summary -----------------------------------------------
cat("Final dataset:\n")
cat("Observations:", nrow(model_data_clean), "\n")
cat("Sites:", length(unique(model_data_clean$site_id)), "\n")
cat("Size category distribution:\n")
print(table(model_data_clean$size_category))
cat("Accuracy distribution:\n")
print(table(model_data_clean$accuracy))

# 6. Bayesian model with Student-t distribution --------------------------
interaction_model_t <- brm(
  bf(log_count ~ centered_year * size_category +
       (centered_year | site_id) +
       (centered_year | accuracy_factor)),
  data = model_data_clean,
  family = student(),
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(2), class = "sd", group = "site_id"),
    prior(student_t(3, 0, 1), class = "sd", group = "accuracy_factor"),
    prior(lkj(4), class = "cor"),
    prior(gamma(2, 0.1), class = "nu")   # degrees of freedom
  ),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)

# Save model object 
#saveRDS(interaction_model_t, "interaction_model_t.rds")

# 7. Model diagnostics ---------------------------------------------------
summary(interaction_model_t)
bayes_R2(interaction_model_t)



# 
# TRACE PLOTS (fuzzy caterpillar plots) for convergence check
# 
# Extract posterior draws as a matrix
posterior_draws <- as_draws_matrix(interaction_model_t)

# Choose the most important parameters to monitor
trace_params <- c(
  "b_Intercept",
  "b_centered_year",
  "b_size_categorysmall",
  "b_size_categorytypical",
  "b_size_categorylarge",
  "b_centered_year:size_categorysmall",
  "b_centered_year:size_categorytypical",
  "b_centered_year:size_categorylarge",
  "sd_site_id__Intercept",
  "sd_site_id__centered_year",
  "sigma",
  "nu"   # degrees of freedom of the Student-t
)

# Create the trace plot
trace_plot <- mcmc_trace(posterior_draws,
                         pars = trace_params,
                         facet_args = list(ncol = 3, scales = "free")) +
  scale_color_manual(values = c("#1b7837", "#7fbf7b", "#ffffbf", "#fc8d59")) +
  labs(title = "Trace plots (fuzzy caterpillar plots)",
       subtitle = "Each colour represents a different MCMC chain. Overlapping, stationary chains indicate good convergence.",
       x = "Iteration",
       y = "Parameter value") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"))

# Display the plot
print(trace_plot)

# Save as PDF (high resolution)
##ggsave("trace_plots.pdf", trace_plot, width = 12, height = 10, dpi = 600)



# Posterior predictive check
y_rep <- posterior_predict(interaction_model_t, ndraws = 100)
ppc_dens_overlay(y = interaction_model_t$data$log_count, yrep = y_rep) +
  labs(title = "Posterior predictive check") +
  theme_bw()+xlim(-20,20)

# Save as PDF (high resolution)
##ggsave("posterior_prediction_plots.pdf", width = 8, height = 6, dpi = 600)



# LOO with moment matching (requires save_pars = TRUE)
loo_result <- loo(interaction_model_t, moment_match = TRUE)
print(loo_result)




# 8. Extract accuracy effects --------------------------------------------
accuracy_intercept <- interaction_model_t %>%
  spread_draws(r_accuracy_factor[accuracy, term]) %>%
  filter(term == "Intercept") %>%
  group_by(accuracy) %>%
  summarise(
    intercept_adj = mean(r_accuracy_factor),
    intercept_lower = quantile(r_accuracy_factor, 0.025),
    intercept_upper = quantile(r_accuracy_factor, 0.975)
  ) %>%
  mutate(accuracy = factor(accuracy, levels = 1:5))

accuracy_slope <- interaction_model_t %>%
  spread_draws(r_accuracy_factor[accuracy, term]) %>%
  filter(term == "centered_year") %>%
  group_by(accuracy) %>%
  summarise(
    slope_adj = mean(r_accuracy_factor),
    slope_lower = quantile(r_accuracy_factor, 0.025),
    slope_upper = quantile(r_accuracy_factor, 0.975)
  ) %>%
  mutate(accuracy = factor(accuracy, levels = 1:5))

accuracy_effects <- accuracy_intercept %>%
  left_join(accuracy_slope, by = "accuracy")

print(accuracy_effects)


ggplot(accuracy_effects,aes(as.factor(accuracy),intercept_adj))+
  geom_errorbar(aes(ymin=intercept_lower,ymax=intercept_upper),width=0.1)+
  geom_point(size=3)+
  theme_bw()

ggplot(accuracy_effects,aes(as.factor(accuracy),slope_adj))+
  geom_errorbar(aes(ymin=slope_lower,ymax=slope_upper),width=0.1)+
  geom_point(size=3)+
  theme_bw()


summary(accuracy_effects$intercept_adj)
summary(accuracy_effects$slope_adj)

# 1. Scale intercept adjustments to [0,1]
intercept_scaled <- accuracy_effects %>%
  select(accuracy, intercept_adj, intercept_lower, intercept_upper) %>%
  mutate(
    intercept_adj_scaled = (intercept_adj - min(intercept_adj)) / (max(intercept_adj) - min(intercept_adj)),
    intercept_lower_scaled = (intercept_lower - min(intercept_adj)) / (max(intercept_adj) - min(intercept_adj)),
    intercept_upper_scaled = (intercept_upper - min(intercept_adj)) / (max(intercept_adj) - min(intercept_adj))
  )

# 2. Scale slope adjustments to [0,1]
slope_scaled <- accuracy_effects %>%
  select(accuracy, slope_adj, slope_lower, slope_upper) %>%
  mutate(
    slope_adj_scaled = (slope_adj - min(slope_adj)) / (max(slope_adj) - min(slope_adj)),
    slope_lower_scaled = (slope_lower - min(slope_adj)) / (max(slope_adj) - min(slope_adj)),
    slope_upper_scaled = (slope_upper - min(slope_adj)) / (max(slope_adj) - min(slope_adj))
  )

# 3. Combine into long format for ggplot
accuracy_combined <- bind_rows(
  intercept_scaled %>%
    select(accuracy, 
           adj = intercept_adj_scaled, 
           lower = intercept_lower_scaled, 
           upper = intercept_upper_scaled) %>%
    mutate(effect = "Intercept"),
  slope_scaled %>%
    select(accuracy, 
           adj = slope_adj_scaled, 
           lower = slope_lower_scaled, 
           upper = slope_upper_scaled) %>%
    mutate(effect = "Slope")
)

# 4. Plot both on the same [0,1] scale
ggplot(accuracy_combined, aes(x = accuracy, y = adj, color = effect)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.3)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c("Intercept" = "#56B4E9", "Slope" = "#E69F00")) +
  labs(
    title = "Accuracy effects: intercept and slope adjustments (separately scaled to [0,1])",
    x = "Accuracy level",
    y = "Normalized adjustment (unitless, within each effect type)",
    color = "Effect"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

##ggsave("accuracy_effects_plots.pdf", width = 8, height = 6, dpi = 600)


# 9. Site-level random effects and correlation ---------------------------
site_effects <- interaction_model_t %>%
  spread_draws(r_site_id[site_id, term]) %>%
  mutate(effect_type = ifelse(term == "Intercept", "intercept", "slope")) %>%
  group_by(site_id, effect_type) %>%
  summarise(
    mean = mean(r_site_id),
    lower = quantile(r_site_id, 0.025),
    upper = quantile(r_site_id, 0.975),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = site_id,
    names_from = effect_type,
    values_from = c(mean, lower, upper),
    names_glue = "{effect_type}_{.value}"
  )




# Add size category and colony size info
colony_info <- model_data_clean %>%
  group_by(site_id, size_category) %>%
  
  
  summarise(
    first_year = min(season_starting, na.rm = TRUE),
    last_year = max(season_starting, na.rm = TRUE),
    first_count = first(penguin_count[season_starting == first_year]),
    first_accuracy = first(accuracy[season_starting == first_year]),
    last_count = first(penguin_count[season_starting == last_year]),
    last_accuracy = first(accuracy[season_starting == last_year]),
    time_interval = last_year - first_year,
    total_change = ifelse(first_count > 0, ((last_count / first_count) - 1) * 100, NA),
    annual_change = total_change / time_interval,
    mean_count = mean(penguin_count, na.rm = TRUE),
    n_years = n_distinct(season_starting),
    lon=mean(longitude_epsg_4326),
    lat=mean(latitude_epsg_4326),
    .groups = 'drop'
  ) %>%
  mutate(total_pop = sum(mean_count),
         weight = mean_count / total_pop)
  

site_effects <- site_effects %>%
  left_join(colony_info, by = "site_id")

head(site_effects)

# Extract correlation between intercept and slope
vc <- VarCorr(interaction_model_t)
vc$accuracy_factor
vc$site_id$cor

ranef_summary <- summary(interaction_model_t)$random$site_id
cor_est <- ranef_summary["cor(Intercept,centered_year)", "Estimate"]
cor_lower <- ranef_summary["cor(Intercept,centered_year)", "l-95% CI"]
cor_upper <- ranef_summary["cor(Intercept,centered_year)", "u-95% CI"]


# 10. Derive percent change over 3 generations (28.2 years) -------------
# For each site, compute total slope = b_centered_year + interaction + r_site_id
population_draws <- interaction_model_t %>%
  spread_draws(b_centered_year,
               `b_centered_year:size_categorysmall`,
               `b_centered_year:size_categorytypical`,
               `b_centered_year:size_categorylarge`,
               r_site_id[site_id, term]) %>%
  filter(term == "centered_year") %>%
  left_join(colony_info, by = "site_id") %>%
  mutate(
    site_slope = case_when(
      size_category == "tiny"     ~ b_centered_year + r_site_id,
      size_category == "small"    ~ b_centered_year + `b_centered_year:size_categorysmall` + r_site_id,
      size_category == "typical"  ~ b_centered_year + `b_centered_year:size_categorytypical` + r_site_id,
      size_category == "large"    ~ b_centered_year + `b_centered_year:size_categorylarge` + r_site_id,
      TRUE ~ b_centered_year + r_site_id
    )
  )

summary(population_draws)


# Population-weighted metrics
population_risk <- population_draws %>%
  group_by(.draw) %>%
  summarise(
    weighted_slope = sum(site_slope * weight, na.rm = TRUE),
    pop_pct_change = (exp(weighted_slope * 3 * GENERATION_LENGTH) - 1) * 100,
    total_increase = sum(mean_count[site_slope > 0] * (exp(site_slope[site_slope > 0] * 3 * GENERATION_LENGTH) - 1), na.rm = TRUE),
    total_decline = sum(abs(mean_count[site_slope < 0] * (exp(site_slope[site_slope < 0] * 3 * GENERATION_LENGTH) - 1)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    net_change = total_increase - total_decline,
    compensation_ratio = total_increase / total_decline
  )

# 11. Summarize population-level risk ------------------------------------
n_draws <- nrow(population_risk)

pop_summary <- population_risk %>%
  summarise(
    median_pop_change = median(pop_pct_change),
    median_net_change= median(net_change),
    lower_pop_change = quantile(pop_pct_change, 0.025),
    upper_pop_change = quantile(pop_pct_change, 0.975),
    
    median_net_change = median(net_change),
    net_change_lower = quantile(net_change, 0.025),
    net_change_upper = quantile(net_change, 0.975),
    
    n_risk_30 = sum(pop_pct_change <= -30),
    n_risk_50 = sum(pop_pct_change <= -50),

    n_compensation = sum(net_change > 0)
  ) %>%
  mutate(
    risk_30pct = n_risk_30 / n_draws,
    risk_30_lower = qbeta(0.025, n_risk_30, n_draws - n_risk_30 + 1),
    risk_30_upper = qbeta(0.975, n_risk_30 + 1, n_draws - n_risk_30),
    risk_50pct = n_risk_50 / n_draws,
    risk_50_lower = qbeta(0.025, n_risk_50, n_draws - n_risk_50 + 1),
    risk_50_upper = qbeta(0.975, n_risk_50 + 1, n_draws - n_risk_50),
    prob_compensation = n_compensation / n_draws,
    comp_lower = qbeta(0.025, n_compensation, n_draws - n_compensation + 1),
    comp_upper = qbeta(0.975, n_compensation + 1, n_draws - n_compensation)
  )

print(pop_summary)

# 12. Create final risk table --------------------------------------------
risk_table <- data.frame(
  Metric = c("Population trend (% change)",
             "Risk of ≥30% decline",
             "Risk of ≥50% decline",
             "Probability of compensation",
             "Net population change (nests)"),
  Estimate = c(
    paste0(round(pop_summary$median_pop_change, 1), "% [",
           round(pop_summary$lower_pop_change, 1), " to ",
           round(pop_summary$upper_pop_change, 1), "]"),
    paste0(round(pop_summary$risk_30pct * 100, 1), "% [",
           round(pop_summary$risk_30_lower * 100, 1), " to ",
           round(pop_summary$risk_30_upper * 100, 1), "]"),
    paste0(round(pop_summary$risk_50pct * 100, 1), "% [",
           round(pop_summary$risk_50_lower * 100, 1), " to ",
           round(pop_summary$risk_50_upper * 100, 1), "]"),
    paste0(round(pop_summary$prob_compensation * 100, 1), "% [",
           round(pop_summary$comp_lower * 100, 1), " to ",
           round(pop_summary$comp_upper * 100, 1), "]"),
    paste0(round(pop_summary$median_net_change, 0), " [",
           round(pop_summary$net_change_lower, 0), " to ",
           round(pop_summary$net_change_upper, 0), "]")
  )
)

print(risk_table)

####------- size category and risk ------------


# Compute size_risk from population_draws
size_risk <- population_draws %>%
  group_by(size_category, .draw) %>%
  summarise(
    category_pct = (exp(mean(site_slope) * 3 * GENERATION_LENGTH) - 1) * 100,
    category_weight = sum(weight),
    .groups = "drop"
  ) %>%
  group_by(size_category) %>%
  summarise(
    median_pct = median(category_pct),
    lower_pct = quantile(category_pct, 0.025),
    upper_pct = quantile(category_pct, 0.975),
    risk_30 = mean(category_pct <= -30),
    risk_50 = mean(category_pct <= -50),
    contribution = mean(category_weight),
    .groups = "drop"
  )



# Create the table
size_risk_table <- size_risk %>%
  mutate(
    size_category = factor(size_category,
                           levels = c("tiny", "small", "typical", "large"),
                           labels = c("Tiny (<80)", "Small (80 to 600)",
                                      "Typical (600 to 9,000)", "Large (>9,000)")),
    `Population share` = paste0(round(contribution * 100, 1), "%"),
    `Median % change [95% CI]` = paste0(
      round(median_pct, 1), "% [",
      round(lower_pct, 1), "% to ",
      round(upper_pct, 1), "%]"
    ),
    `Risk of ≥30% decline` = paste0(round(risk_30 * 100, 1), "%"),
    `Risk of ≥50% decline` = paste0(round(risk_50 * 100, 1), "%")
  ) %>%
  select(size_category, `Population share`, `Median % change [95% CI]`,
         `Risk of ≥30% decline`, `Risk of ≥50% decline`)

# Print the table
print(size_risk_table)

# Save as CSV
#write.csv(size_risk_table, "size_category_risk_table.csv", row.names = FALSE)

#write.csv(risk_table,"risk_table.csv")

# 13. Generate final figures  ----------------------------------

head(site_effects)
# Figure: Site-level random effects correlation
p_corr <- site_effects %>%
  ggplot(aes(x = intercept_mean, y = slope_mean, color = size_category)) +
  geom_errorbar(aes(ymin = slope_lower, ymax = slope_upper), alpha = 0.3, width = 0) +
  geom_errorbarh(aes(xmin = intercept_lower, xmax = intercept_upper), alpha = 0.3, height = 0) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", alpha = 0.2) +
  scale_color_manual(values = c("tiny" = "steelblue", "small" = "green3",
                                "typical" = "orange3", "large" = "red3")) +
  labs(title = "Site-level random effects",
       subtitle = paste0("Correlation = ", round(cor_est, 2),
                         " [", round(cor_lower, 2), "-", round(cor_upper, 2), "]"),
       x = "Intercept adjustment (baseline size deviation)",
       y = "Slope adjustment (trend deviation)",
       color = "Size category") +
  theme_bw() +
  theme(legend.position = "bottom")

p_corr

##ggsave("site_correlation.pdf", p_corr, width = 8, height = 6)

med_net_change<-median(population_risk$net_change)

# Figure: Net change density (compensation)
p_comp <- population_risk %>%
  ggplot(aes(x = net_change)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = med_net_change, linetype = "dashed", color = "red3") +
  labs(x = "Net change in number of nests",
       y = "Density",
       title = "Net population change over 3 generations") +
  scale_x_continuous(labels = scales::comma,limits =c(-500000,200000)) +
  
  theme_bw()+theme(axis.text.y=element_blank())

p_comp

#ggsave("compensation_plots.pdf", width = 6, height = 4, dpi = 600)

# 14. Save final tables --------------------------------------------------
##write.csv(accuracy_effects, "accuracy_effects.csv", row.names = FALSE)
##write.csv(risk_table, "risk_table_final.csv", row.names = FALSE)
##write.csv(pop_summary, "population_risk_summary.csv", row.names = FALSE)


# ---------- site effects geographically ---------------



library(sf)
library(dplyr)
library(terra)

land<-st_as_sf(vect("C:/GIS/add_coastline_medium_res_polygon_v7_8.shp/add_coastline_medium_res_polygon_v7_8.shp")
)

land<-st_transform(land,4326)
plot(land)


d1mpa<-st_as_sf(vect("C:/CCAMLR 2024/SC/StephanVersion/HS2024_D1MPA/HS_2024_D1MPA_6932.shp"))%>%
  st_transform(4326)
plot(d1mpa)

table(d1mpa$cat_sc08)

d1mpa$categ[d1mpa$cat_sc08!="GPZ" & d1mpa$cat_sc08!="TRZ" ]<-"SPZ"

table(d1mpa$categ)

d1mpa$categ[d1mpa$categ=="TRZ"]<-"GPZ"

table(d1mpa$categ)
head(site_effects)


site_effects<-site_effects%>%
  mutate(magnitude=abs(slope_mean),
         trend=ifelse(slope_mean<0,"-","+"))

ggplot() +
 # geom_sf(data = d1mpa, aes(fill = categ), alpha = 0.5) +
  geom_sf(data = land, fill = "grey80", colour = "grey40") +
  geom_point(data = site_effects, aes(lon, lat, size = mean_count,
                                      
                            fill = (slope_mean), shape = trend),
             alpha=0.75) +
  scale_fill_gradient2(high = "blue", mid = "grey90", low= "red3",
                         midpoint = 0, name = "Slope") +

  scale_size_continuous(range = c(1, 10), name = "Colony size") +
  scale_shape_manual(values = c("+" = 21, "-" = 25)) +
  coord_sf(xlim = c(-70, -40), ylim = c(-68, -60)) +
  theme_bw() +
  
  theme(legend.position = "right")


# classify sites based on regions

site_effects$location[site_effects$lon>(-50)]<-"SOI"

site_effects$location[site_effects$lon<(-50) &
                        site_effects$lon>(-56.75) ]<-"EI"

site_effects$location[is.na(site_effects$location) &
                        site_effects$lat>(-63.3) ]<-"SSI"

site_effects$location[is.na(site_effects$location) ]<-"AP"

table(site_effects$location)

site_effects$location<-factor(site_effects$location,
                              levels=c("AP","SSI","EI","SOI"))

# check if it is right

ggplot() +
  # geom_sf(data = d1mpa, aes(fill = categ), alpha = 0.5) +
  geom_sf(data = land, fill = "grey80", colour = "grey40") +
  geom_point(data = site_effects, aes(lon, lat, 
                                      
                                      colour = location, shape = location),
             alpha=0.75)+
  coord_sf(xlim = c(-70, -40), ylim = c(-68, -60))


ggplot(site_effects,aes(location,magnitude,fill=trend))+geom_boxplot()+
  ylim(0,0.12)+
  
  
ggplot(site_effects,aes(location,mean_count,fill=trend))+geom_boxplot()+
  scale_y_log10(labels=scales::comma_format())
 
ggplot(site_effects,aes(location,n_years,fill=trend))+geom_boxplot()+
  scale_y_log10(labels=scales::comma_format())



ggplot(site_effects,aes(mean_count,magnitude,colour=trend,shape=trend,
                        linetype=trend))+
  geom_smooth(method="lm",se=F)+
  geom_point(aes(size=n_years))+
  scale_x_log10(labels=scales::comma_format())+
  facet_wrap(location~.,scales="free")+
  theme_bw()+
  scale_colour_manual(values=c("red3","steelblue"),name="Trend") +
  scale_linetype_manual(values=c("solid","dashed"),name="Trend") +
  scale_size_continuous(range = c(1, 10), name = "Counts") +
  scale_shape_manual(values = c("+" = 21, "-" = 25),name="Trend") +
  xlab("Mean colony size (nests)")+
  ylab("|Slope|")
  



#ggsave("colonies_trends.pdf", width = 9, height = 7, dpi = 600)

# End of script