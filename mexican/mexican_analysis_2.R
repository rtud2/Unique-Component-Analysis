setwd("mexican")
temp2 <- fread("merged_all_gender.geno", header = F) # 2297 SNPs x 207 Sample ID
temp2_names <- fread("merged_all_gender_names", header = F) #Sample ID

# transforming data into a usable data.table
geno_to_R = function(data_table){
  data_table <- data_table[, tstrsplit(V1, "", fixed=TRUE)][,lapply(.SD, as.numeric)]
}
temp2 <- geno_to_R(temp2)

#creating a lookup table to figure out which ethnic group the samples belong to
mex_group <- c("GUA","GUE","SON", "VER", "YUC", "ZAC")
mex_grp_names <- c(sapply(mex_group, function(xx) paste0(xx ,"_", 1:50)), paste0("ZAP_", 1:30), paste0("CEU", 1:90))

setnames(temp2, mex_grp_names)



t_temp <- data.table(t(temp2))

# which SNPs have missing values in them
imputed_cols <- t_temp[, which(lapply(.SD, function(zz) sum(zz == 9)) > 0)]

#index at those SNPs that have missing values
list_missing_loc <- t_temp[,.SD, .SDcols = imputed_cols][, .(lapply(.SD, function(zz) which(zz == 9)))][[1]]

for(jj in seq_along(imputed_cols)){
  idx_replace <- list_missing_loc[[jj]]
  which_col <-  names(imputed_cols[jj])
  
  #getting the mean for the particular SNP without the missing values
  all_snp_vals <- t_temp[, unlist(.SD), .SDcols = which_col]
  
  #replacing the index 
  t_temp[idx_replace, (which_col):= mean(all_snp_vals[which(all_snp_vals!=9)]) ]
}

imputed_mean_dat <- data.table(t(t_temp))

setnames(imputed_mean_dat, mex_grp_names)

not_ZAP <- grep(mex_grp_names, pattern = "ZAP|CEU", invert = T,   value = T)
ZAP <- grep(mex_grp_names, pattern = "ZAP|CEU",   value = T)
length(not_ZAP)

#before imputation: 
scaled_temp <- temp[, scale(t(.SD)), .SDcols = not_ZAP]
scaled_temp <- scaled_temp[, which(!is.nan(colSums(scaled_temp)))] #remove columns with NA values

centered_temp <- temp[, scale(t(.SD), scale = F), .SDcols = not_ZAP]
centered_temp <- centered_temp[, which(!is.nan(colSums(scaled_temp)))] #remove same columns as scaled

scaled_rotated_dat <- scaled_temp %*% svd(scaled_temp, nv = 2)$v
centered_rotated_dat <- centered_temp %*% svd(centered_temp, nv = 2)$v

plot_dat <- rbind(data.table(centered_rotated_dat, substr(not_ZAP, start = 1, stop =3), "Unimputed","Centered Only"),
                  data.table(scaled_rotated_dat, substr(not_ZAP, start = 1, stop =3), "Unimputed","Scaled"))
setnames(plot_dat, c("PC1", "PC2", "State", "Imputation","Scale"))


#after imputation: distribution imputation
#scaled_imp_dist <- imputed_dat[, scale(t(.SD)), .SDcols = not_ZAP]
#scaled_imp_dist <- scaled_imp_dist[, which(!is.nan(colSums(scaled_imp_dist)))] #remove columns with NA values

#centered_imp_dist <- imputed_dat[, scale(t(.SD), scale = F), .SDcols = not_ZAP]
#centered_imp_dist <- centered_imp_dist[, which(!is.nan(colSums(scaled_imp_dist)))] #remove same columns as scaled

#after imputation: mean imputation
 scaled_imp_mean <- imputed_mean_dat[, scale(t(.SD)), .SDcols = not_ZAP]
 scaled_imp_mean <- scaled_imp_mean[, which(!is.nan(colSums(scaled_imp_mean)))] #remove columns with NA values
# 
centered_imp_mean <- imputed_mean_dat[, scale(t(.SD), scale = F), .SDcols = not_ZAP]
centered_imp_mean <- centered_imp_mean[, which(!is.nan(colSums(scaled_imp_mean)))] #remove same columns as scaled
# 
# imp_dist_scaled_rotated_dat <- scaled_imp_dist %*% svd(crossprod(scaled_imp_dist), nv = 2)$v
# imp_dist_centered_rotated_dat <- centered_imp_dist %*% svd(centered_imp_dist, nv = 2)$v

imp_mean_scaled_rotated_dat <- scaled_imp_mean %*% svd(scaled_imp_mean, nv = 2)$v
imp_mean_centered_rotated_dat <- centered_imp_mean %*% svd(centered_imp_mean, nv = 2)$v

plot_imp_dat <- rbind(#data.table(imp_dist_scaled_rotated_dat, substr(not_ZAP, start = 1, stop =3), "Dist.","Scaled"),
                      #data.table(imp_dist_centered_rotated_dat, substr(not_ZAP, start = 1, stop =3), "Dist.","Centered Only"),
                      data.table(imp_mean_scaled_rotated_dat, substr(not_ZAP, start = 1, stop =3), "Mean","Scaled"),
                      data.table(imp_mean_centered_rotated_dat, substr(not_ZAP, start = 1, stop =3), "Mean","Centered Only"))
setnames(plot_imp_dat, c("PC1", "PC2", "State", "Imputation","Scale"))


#plot both side-by-side
comb_plot_dat <- rbind(plot_dat, plot_imp_dat)

pre_post_imputation_pca <- ggplot(data = comb_plot_dat)+
  geom_jitter(aes(x = PC1, y = PC2, color = State))+
  facet_grid(Scale~Imputation)+
  labs(title = "Principal Components of Mexican Heritage")

ggsave(plot = pre_post_imputation_pca, filename = "pre_post_imputation_pca.png",height = 8.5, width = 11, units = "in")