library(sf)            # spatial data handling
library(spdep)         # spatial weights, tests
library(spatialreg)    # errorsarlm, lagsarlm
# spatialprobit provides spatial probit/EM; optional but recommended
# If not installed, the script falls back to non-spatial glm + spaMM suggestion
suppressWarnings(suppressMessages(library(spatialprobit)))
library(ggplot2)
library(tmap)
library(dplyr)
library(texreg)
library(broom)




##################################################################
# SHAPEFILE CREATION AND VISUALIZATION
mydata <- read.csv("cocoa.csv")
str(mydata)

mydata_sf <- st_as_sf(mydata, coords = c("longitude", "latitude"))
str(mydata_sf)


#ASSIGN WGS84 PROJECTION
st_crs(mydata_sf) <- 4326
#CHECK
st_crs(mydata_sf)

mydata_sf
head(mydata_sf)
sf::st_geometry_type(mydata_sf)

# 
# 
# #BASEMAP
# library(ggmap)
# citation("ggmap")
# #REGISTER  STADIA MAP FOR API
# register_stadiamaps(key = "aa4a6990-826c-4e2c-8525-862b0318eb09")
# 
# mydata_sf %>% st_bbox() %>% as.vector() %>% 
#  ggmap::get_stadiamap(zoom = 11, messaging = FALSE) -> baseMap;
# 
# p <- ggmap(baseMap) + 
#   geom_point(aes(x=X, y=Y), 
#              data = mydata_sf %>% st_coordinates() %>% tibble :: as_tibble(),
#              color = "brown",
#              size = 1,
#              alpha = 0.5
#              )
# ggsave(
#   filename = "high_res_map.png",
#   plot = p,
#   width = 10,        # in inches
#   height = 8,        # in inches
#   dpi = 400          # typical for print-quality
# )


#To convert sf to sp just use, sf::st_as_sf
#Save the sf
st_write(mydata_sf, "cocoa.shp", delete_layer = TRUE)
#############################################################################

#READ NEW DATA
cocoa_sf <- st_read("cocoa.shp", quiet=TRUE)
head(cocoa_sf)
par(mar=c(1,1,1,1))
plot(cocoa_sf)

dev.off()
#Plot of yield per hectare
plot(cocoa_sf$yld_pr_h)
hist(cocoa_sf$yld_pr_h, axes = TRUE, ann = TRUE)

#BETTER VISUALIZE/DISTINGUISH YIELDS WITH LOW VALUES
ggplot(cocoa_sf) +
  geom_sf(aes(color = yld_pr_h)) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    title = "Cocoa Yield Distribution",
    subtitle = "Yield values across sampled farm locations",
    color = "Yield (kg/ha)"
  )


#RECHECK THE GEOMETRY AGAIN TO CONFIRM
sf::st_geometry_type(cocoa_sf)


#SPATIAL AUTOCORRELATION (GLOBAL AND LISA)
head(cocoa_sf)
str(cocoa_sf)
names(cocoa_sf)


library(magrittr)

# 1. Extract coordinates from your sf object
coords <- st_coordinates(cocoa_sf)

# 2. Find k nearest neighbors (e.g., k = 4)
knn_obj <- knearneigh(coords, k = 7)

# 3. Convert to neighbors list
nb_obj <- knn2nb(knn_obj)

# 4. Convert to weights list
cocoaNbList <- nb2listw(nb_obj, zero.policy = TRUE)

#Check how neighborhood works
summary(nb_obj)
n.comp.nb(nb_obj)   # number of connected components

# Global Moran's I test
options(scipen = 5)
cocoaNbList %>% 
  spdep::moran.test(cocoa_sf$yld_pr_h, ., zero.policy = TRUE)

# Global Moran's I plot
spdep::moran.plot(cocoa_sf$yld_pr_h,
                  cocoaNbList,
                  zero.policy = TRUE,
                  xlab = 'Cocoa yield per hectare',
                  ylab = 'Lagged cocoa yield per hectare (of Neighbors)',
                  pch = 20)

#Local Moran's I plot
lisarslt <- spdep::localmoran(cocoa_sf$yld_pr_h,
                              cocoaNbList,
                              zero.policy = TRUE,
                              na.action = na.omit
)
dim(lisarslt); dim(cocoa_sf);
head(lisarslt)

#Derive the cluster/outlier types (CoType in ArcGIS) for each spatial feature in the data
significanceLevel <- 0.05  # 95% confidence
meanVal <- mean(cocoa_sf$yld_pr_h)

#For Debugging
#lisarslt <- lisarslt %>% dplyr::select(-coType)
#names(lisarslt)

lisarslt <- lisarslt %>%
  tibble::as_tibble() %>%
  magrittr::set_colnames(c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "Pr")) %>%   # ✅ renamed here
  dplyr::mutate(
    coType = dplyr::case_when(
      Pr > significanceLevel ~ "Insignificant",
      Pr <= significanceLevel & Ii >= 0 & cocoa_sf$yld_pr_h >= meanVal ~ "HH",
      Pr <= significanceLevel & Ii >= 0 & cocoa_sf$yld_pr_h <  meanVal ~ "LL",
      Pr <= significanceLevel & Ii <  0 & cocoa_sf$yld_pr_h >= meanVal ~ "HL",
      Pr <= significanceLevel & Ii <  0 & cocoa_sf$yld_pr_h <  meanVal ~ "LH",
      TRUE ~ "Undefined"
    )
  )

# Attach classification
cocoa_sf$lis_cluster <- lisarslt$coType

# Plot map
ggplot(cocoa_sf) +
  geom_sf(aes(color = lis_cluster)) +
  scale_color_manual(
    values = c(
      "HH" = "red",
      "LL" = "blue",
      "HL" = "orange",
      "LH" = "green",
      "Insignificant" = "grey80"
    )
  ) +
  theme_minimal() +
  labs(
    title = "Local Moran's I Cluster Map",
    subtitle = "Cocoa Yield Spatial Clusters",
    color = "Cluster Type"
  )




#LOGISTIC REGRESSION TO TEST PRODUCTIVITY
cocoa_sf <- cocoa_sf %>%
  mutate(yield_bin = ifelse(yld_pr_h > 0, 1L, 0L))

#Number of zeros(no yields), number of ones(yields)
table(cocoa_sf$yield_bin)

pa = yield_bin ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
  fdp_prv + map_frm + gender + hshld_s + nb_frms + yld_pr_t + pds_pr_ + tr_dnst



#Convert the weight list to matrix as part of spatial probit requirement
cocoaNbList_mat <- nb2mat(nb_obj, style = "W", zero.policy = TRUE)
cocoaNbList_mat <- Matrix(cocoaNbList_mat, sparse = TRUE)
diag(cocoaNbList_mat) <- 0


#Model1: SAR Spatial Probit model using MCMC
sar_pm <- spatialprobit::sarprobit(
  formula = pa,
  data    = cocoa_sf,
  W       = cocoaNbList_mat,
  ndraw   = 5000,        # start with 1000 draws; increase to 5000 for final model
  showProgress = TRUE
)

summary(sar_pm)












#PART B (NOW MODELLING YIELD WITH FARMS WITH POSITIVE YIELDS)
#Extract only yields > 0
# Define predictor names (vector of variable names)
predictors <- c("farm_sz", "srfc_mp", "nb_f_2122", "rcvd_rm", 
                "c_21_22", "d_frmr_", "fdp_prv", "map_frm", 
                "gender", "hshld_s", "nb_frms", "yld_pr_t", 
                "pds_pr_", "tr_dnst")

# Identify rows with positive yield and complete predictor data
head(cocoa_sf)
# Drop geometry temporarily
cocoa_df <- st_drop_geometry(cocoa_sf)
head(cocoa_df)

# Now safely subset
pos_yield <- which(cocoa_df$yld_pr_h > 0 & 
                     complete.cases(cocoa_df[, predictors]))

# Subset the original spatial data
cocoa_pos <- cocoa_sf[pos_yield, ]
head(cocoa_pos)

# Check
nrow(cocoa_pos)

#SPATIAL REGRESSION
names(cocoa_pos)

#Model1
library(car)

ols1 = lm(yld_pr_h ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
            fdp_prv + map_frm + gender + hshld_s + nb_frms + yld_pr_t + pds_pr_ + tr_dnst,
          data = cocoa_pos)


bc <- powerTransform(ols1)
summary(bc)

ols1_trans <- update(ols1, formula = I(yld_pr_h ^ 0.39) ~ .)
summary(ols1_trans)


par(mfrow = c(1, 2))
plot(ols1, which = 1, main = "Original Model")
plot(ols1_trans, which = 1, main = "Transformed Model (λ = 0.39)")

# Cook’s Distance plot
plot(ols1_trans, which = 4, cook.levels = c(4/length(ols1_trans$fitted.values)),
     main = "Cook's Distance (Influence of Observations)")


influencePlot(ols1_trans, 
              id.method = "identify", 
              main = "Influence Plot: Studentized Residuals vs Leverage",
              sub = "Circle size ∝ Cook’s Distance")

plot(hatvalues(ols1_trans), cooks.distance(ols1_trans),
     xlab = "Leverage", ylab = "Cook's Distance",
     main = "Leverage vs. Cook's Distance")
abline(h = 4/length(ols1_trans$fitted.values), col = "red", lty = 2)
abline(v = 2*mean(hatvalues(ols1_trans)), col = "blue", lty = 2)


# Define Cook's distance threshold
n <- nrow(cocoa_pos)
cutoff <- 4 / n

# Find influential points
influential_points <- which(cooks.distance(ols1_trans) > cutoff)

# Convert to regular data frame
infl_df <- as.data.frame(influence.measures(ols1_trans)$infmat)
colnames(infl_df)


# View influence stats for those points
infl_df[influential_points, c("cook.d", "hat", "dffit", "cov.r")]
cocoa_pos[influential_points, ]

# Drop the row and reset row numbering
cocoa_clean <- cocoa_pos[-1269, ]
rownames(cocoa_clean) <- NULL

ols_no1269 <- lm(formula(ols1_trans), data = cocoa_clean)
any(rownames(ols_no1269$model) == "1269")
which.max(cooks.distance(ols_no1269))

# Compare coefficients
summary(ols1_trans)
summary(ols_no1269)

plot(ols_no1269)
res <- residuals(ols_no1269)
hist(res)

library(lmtest)
bptest(ols1_trans)
#BP TEST IS SIGNFICANT, FURTHER NECESSITATING THE NEED FOR SPATIAL ERROR MODEL OR MORAN'S I RESIDUAL CHECK ON THE TRANSFORMED OLS


#Moran's residual test for spatial autocorrelation
#Redfine the neighborhood list
coords <- st_coordinates(cocoa_pos)

# 2. Find k nearest neighbors (e.g., k = 4)
knn_obj <- knearneigh(coords, k = 7)

# 3. Convert to neighbors list
nb_obj <- knn2nb(knn_obj)

# 4. Convert to weights list
yieldlist <- nb2listw(nb_obj, zero.policy = TRUE)


# --- Moran's I Test on OLS Residuals ---
ols_residuals <- residuals(ols_no1269)
cat("\n=== Moran's I Test on OLS Residuals ===\n")
moran_test_ols <- moran.test(ols_residuals, yieldlist)
print(moran_test_ols)

# --- Rao's Score Tests for Spatial Dependence (recommended replacement for lm.LMtests) ---
cat("\n=== Rao's Score Tests (lm.RStests) ===\n")
rs_tests <- lm.RStests(ols1_trans, yieldlist)
print(rs_tests)


#running SER model since moran's residual test is significant
serr <- spatialreg::errorsarlm(I(yld_pr_h^0.39) ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
                                 fdp_prv  + map_frm + gender + hshld_s + nb_frms + yld_pr_t + pds_pr_ + tr_dnst,
                               data = cocoa_pos,
                               listw = yieldlist,
                               zero.policy = TRUE,
                               na.action = na.omit);
summary(serr)
#Deriving the residuals from the regression, while making sure there's no missing values
cat("\n=== Moran’s I on SEM Residuals ===\n")
moran_sem <- moran.test(residuals(serr), yieldlist)
print(moran_sem)

#SAR Model
sarr <- spatialreg::lagsarlm(I(yld_pr_h^0.39) ~ farm_sz + srfc_mp + nb_f_2122 + rcvd_rm +  c_21_22 + d_frmr_ + 
                               fdp_prv  + map_frm + gender + hshld_s + nb_frms + yld_pr_t + pds_pr_ + tr_dnst,
                             data = cocoa_pos,
                             listw = yieldlist,
                             zero.policy = TRUE,
                             na.action = na.omit);
summary(sarr)

#Deriving the residuals from the regression, while making sure there's no missing values

cat("\n=== Moran’s I on SAR Residuals ===\n")
moran_sem <- moran.test(residuals(sarr), yieldlist)
print(moran_sem)

#MODEL DIAGNOSTICS
# Get fitted and residuals from the spatial error model
cocoa_pos$predicted <- fitted(serr)
cocoa_pos$residuals <- residuals(serr)

# Back-transform predicted and residual values
cocoa_pos$predicted_bt <- cocoa_pos$predicted^(1 / 0.39)
cocoa_pos$residual_bt  <- cocoa_pos$yld_pr_h - cocoa_pos$predicted_bt

#ACTUAL VS PREDICTED
# --- Actual vs Predicted Plot (Back-Transformed SER Model) ---

# Keep only complete pairs
valid <- complete.cases(cocoa_pos$yld_pr_h, cocoa_pos$predicted_bt)
obs <- cocoa_pos$yld_pr_h[valid]
pred <- cocoa_pos$predicted_bt[valid]

# Compute model fit statistics
rmse_bt <- sqrt(mean((obs - pred)^2))
r2_bt <- cor(obs, pred)^2  # correlation-based R²

# Plot
plot(obs, pred,
     xlab = "Observed Yield (kg/ha)",
     ylab = "Predicted Yield (kg/ha)",
     main = "Actual vs Predicted Cocoa Yield (Back-Transformed)",
     pch = 19, col = "darkblue",
     cex = 0.9)

# Add 1:1 line (perfect fit)
abline(0, 1, col = "red", lwd = 2)

# Add fitted regression line
abline(lm(pred ~ obs), col = "green", lwd = 2, lty = 2)

# Add legend
legend("topleft",
       legend = c("1:1 Line", "Fitted Line"),
       col = c("red", "green"),
       lwd = 2,
       lty = c(1, 2),
       bty = "n")

# Add annotation with model metrics
text(x = max(obs, na.rm = TRUE) * 0.55,
     y = max(pred, na.rm = TRUE) * 0.45,
     labels = paste("R² =", round(r2_bt, 3),
                    "\nRMSE =", round(rmse_bt, 2), "kg/ha"),
     adj = 0, cex = 0.9, font = 2)

library(tmap)
tmap_mode("plot")

tm_shape(cocoa_pos) +
  tm_dots(col = "residual_bt",
          palette = "-RdBu",
          midpoint = 0,
          size = 0.2,
          style = "quantile",
          title = "Residuals (kg/ha)") +
  tm_layout(main.title = "Spatial Distribution of Model Residuals (Back-Transformed)",
            main.title.size = 1.1,
            legend.outside = TRUE,
            legend.outside.position = "right") +
  tm_compass(type = "8star", position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom"))

#Model prediction error (Absolute)
cocoa_pos$abs_error_bt <- abs(cocoa_pos$residual_bt)
tm_shape(cocoa_pos) +
  tm_dots(col = "abs_error_bt",
          palette = "YlOrRd",
          size = 0.2,
          style = "quantile",
          title = "Absolute Error (kg/ha)") +
  tm_layout(main.title = "Spatial Distribution of Prediction Errors (Back-Transformed)",
            main.title.size = 1.1,
            legend.outside = TRUE,
            legend.outside.position = "right") +
  tm_compass(type = "8star", position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom"))


# === Model Accuracy Metrics (Back-Transformed Scale) ===
# RMSE and MAE
rmse_bt <- sqrt(mean(cocoa_pos$residual_bt^2, na.rm = TRUE))
mae_bt  <- mean(abs(cocoa_pos$residual_bt), na.rm = TRUE)

# R-squared on back-transformed scale
sst_bt <- sum((cocoa_pos$yld_pr_h - mean(cocoa_pos$yld_pr_h, na.rm = TRUE))^2, na.rm = TRUE)
sse_bt <- sum(cocoa_pos$residual_bt^2, na.rm = TRUE)
r2_bt  <- 1 - (sse_bt / sst_bt)

# Print results nicely
cat("Model accuracy (Back-transformed scale):\n")
cat("RMSE:", round(rmse_bt, 2), "kg/ha\n")
cat("MAE: ", round(mae_bt, 2), "kg/ha\n")
cat("R²:  ", round(r2_bt, 3), "\n")

mean(cocoa_pos$yld_pr_h)

