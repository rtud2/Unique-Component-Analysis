rm(list = ls())
setwd("~/Desktop/Residual-Dimension-Reduction")
#install.packages("uca_0.13.tar.gz", type="source", repos = NULL) 
libraries <- c("data.table", "MASS", "ggplot2", "splines", "gridExtra", "uca", "imager", "future", "RSpectra", "microbenchmark", "Rfast", "rgl", "car")
invisible(lapply(libraries, library, character.only = T))
#plan(tweak(multicore, workers = 8L))
plan(multisession)
## Helper functions 
source("functions/helper_functions.R")
kdef <- "faces/KDEF_and_AKDEF/KDEF/"
kdef_A <- "faces/KDEF_and_AKDEF/A_kdef_files/"

# Practice Emotions
sus_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz, "/", zz, "SUS.JPG")))
sas_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz, "/", zz, "SAS.JPG")))
nes_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz, "/", zz, "NES.JPG")))
dis_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz, "/", zz, "DIS.JPG")))
ans_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz, "/", zz, "ANS.JPG")))
has_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz, "/", zz, "HAS.JPG")))
afs_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz, "/", zz, "AFS.JPG")))

practice_list <- list(sus_pract_files, has_pract_files, sas_pract_files, nes_pract_files, dis_pract_files, ans_pract_files, afs_pract_files)
practice_female_list <- lapply(practice_list, function(zz){grep(pattern = "AF[0-9][0-9]", x = zz, value = T)})
female_emotion_list %<-% lapply(practice_female_list, function(x) readImage(x, 50, 68))
female_emotion_list <- lapply(female_emotion_list, transpose)

# Final Emotions
sus_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz, "/", zz, "SUS.JPG")))
sas_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz, "/", zz, "SAS.JPG")))
nes_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz, "/", zz, "NES.JPG")))
dis_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz, "/", zz, "DIS.JPG")))
ans_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz, "/", zz, "ANS.JPG")))
has_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz, "/", zz, "HAS.JPG")))
afs_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz, "/", zz, "AFS.JPG")))

final_list <- list(sus_final_files, has_final_files, sas_final_files, nes_final_files, dis_final_files, ans_final_files, afs_final_files)
final_female_list <- lapply(final_list, function(zz){grep(pattern = "BF[0-9][0-9]", x = zz, value = T)})
final_female_emotion_list %<-% lapply(final_female_list, function(x) readImage(x, 50, 68))
final_female_emotion_list <- lapply(final_female_emotion_list, transpose)

#### Every Combination of One-Background-Removed ####

all_final <- scale(do.call(rbind, final_female_emotion_list))
kdef_emotions <- c("Surprise", "Happy", "Sad", "Neutral", "Disgust", "Angry", "Afraid")
final_labels <- rep(kdef_emotions, each = 35)
all_pca <- eigs_sym(cov(all_final), k = 5, "LA")
res_list <- list()

#ffor (emo in seq_along(female_emotion_list)) {
for (emo in seq_along(kdef_emotions)) {
    emo_list <- female_emotion_list[-emo]
    emo_list_stack <- as.matrix(scale(do.call(rbind, emo_list)))

    tmp_split <- uca(all_final, lapply(emo_list, scale), nv = 5, method = "data")
    tmp_stack <- uca(all_final, emo_list_stack, nv = 5, method = "data")

    method_labs <- c("PCA", "Pooled", "Split")
    method_obj <- list(all_pca, tmp_stack, tmp_split)
    tmp_projected <- lapply(seq_along(method_obj), function(iter) {
        data.table(all_final %*% method_obj[[iter]]$vectors, final_labels, method_labs[iter])
    })
    tmp_projected <- do.call(rbind, tmp_projected)
    setnames(tmp_projected, c(paste0("V",1:5), "Emotion", "Method"))

    all_emotion_faces <- do.call(rbind, lapply(lapply(method_obj, "[[", "vectors"),
                                function(vectors){toEigenfaces(vectors, 50, 68, "")}))

    all_emotion_faces$alpha <- rep(method_labs, each = 17E3)
    res_list[[emo]] <- list(tmp_projected, all_emotion_faces)
    message(paste0("finished UCA for ", kdef_emotions[emo]))
    saveRDS(res_list[[emo]], paste0("example/Faces/One_Background_Removed/RDS/", kdef_emotions[emo], "removed.rds"))
}

# Saving the data so i dont have to re-run this everytime dave asks me to tweak a figure
# res_list <- lapply(seq_along(kdef_emotions), function(emo){
#                        readRDS(paste0("example/Faces/One_Background_Removed/RDS/", kdef_emotions[emo], "removed.rds"))
#     })

csrsc1 <- c( "#ff852b", "#9f009a", "#90d78e", "#ff5ea1", "#885200", "#90afff", "#7a4527" )
#Plotting the images
for (emo in seq_along(res_list)) {
    tmp_projected <- res_list[[emo]][[1]]
    all_emotion_faces <- res_list[[emo]][[2]]

    tmp_projected_plot <- ggplot(data = tmp_projected) +
    geom_point(aes(x = V1, y = V2, color = Emotion)) +
    facet_wrap(~ Method, scale = "free", nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18)) +
    guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual(values = csrsc1) +
    labs(x = "Component 1", y = "Component 2")
    
    ggsave(paste0("example/Faces/One_Background_Removed/", kdef_emotions[emo], "_projected_1_2_removed.png"),
       tmp_projected_plot,
       width = 11, height = 8, units = "in")

    tmp_projected_plot2 <- ggplot(data = tmp_projected) +
    geom_point(aes(x = V2, y = V3, color = Emotion)) +
    facet_wrap(~ Method, scale = "free", nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18)) +
    guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual(values = csrsc1) +
    labs(x = "Component 2", y = "Component 3")

    ggsave(paste0("example/Faces/One_Background_Removed/all_", kdef_emotions[emo], "_removed.png"),
       tmp_eigenfaces,
       width = 8, height = 11, units = "in")

    message(paste0("finish ", kdef_emotions[emo]))
}

# they're not technically eigenfaces anymore since its the mean
final_mean_faces %<-% lapply(final_female_emotion_list, colMeans)
final_mean_faces_plot <- (toEigenfaces(t(do.call(rbind, final_mean_faces)), 50,68,""))
final_mean_faces_plot[, Emotions := kdef_emotions[Eigenface]]


ggplot(final_mean_faces_plot) +
    geom_raster(aes(x, y, fill = value)) +
    facet_grid(~Emotions) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0),trans = scales::reverse_trans()) +
    scale_fill_gradient(low = "black", high = "white") +
    theme(strip.text = element_text(size = 20),
          axis.title = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank(),
          axis.text =  element_blank())

ggsave("example/Faces/One_Background_Removed/mean_emotions.png",
       width = 8, height = 6, units = "in", dpi = 300)
