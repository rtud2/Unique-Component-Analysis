rm(list = ls())
Sys.setenv(MKL_DEBUG_CPU_TYPE=5)#, MKL_VERBOSE=1) #intel MKL settings. latter to hack amd processors
setwd("/mnt/d/Residual-Dimension-Reduction")
#install.packages("uca_0.13.zip", repos = NULL, type="source")
#install.packages("uca_0.13.tar.gz",type="source", repos = NULL) 
libraries <- c("data.table", "MASS", "ggplot2", "splines","gridExtra", "uca", "imager", "future", "RSpectra", "microbenchmark", "Rfast")

lapply(libraries, library, character.only = T)
plan(multiprocess)

## Helper functions 
source("functions/helper_functions.R")
kdef <- "faces/KDEF_and_AKDEF/KDEF/"
kdef_A <- "faces/KDEF_and_AKDEF/A_kdef_files/"
sus_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"SUS.JPG")))
dis_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"DIS.JPG")))
has_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"HAS.JPG")))
ans_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"ans.JPG")))

sus_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"SUS.JPG")))
dis_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"DIS.JPG")))
afs_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"AFS.JPG")))

#break down into male/female
sus_m_final_files <- grep(pattern = "BM", x = sus_final_files, value = T)
sus_f_final_files <- grep(pattern = "BF", x = sus_final_files, value = T)
dis_m_final_files <- grep(pattern = "BM", x = dis_final_files, value = T)
dis_f_final_files <- grep(pattern = "BF", x = dis_final_files, value = T)
has_m_final_files <- grep(pattern = "BM", x = has_final_files, value = T)
has_f_final_files <- grep(pattern = "BF", x = has_final_files, value = T)
ans_m_final_files <- grep(pattern = "BM", x = ans_final_files, value = T)
ans_f_final_files <- grep(pattern = "BF", x = ans_final_files, value = T)

#background files
sus_f_pract_files <- grep(pattern = "AF", x = sus_pract_files, value = T)
sus_m_pract_files <- grep(pattern = "AM", x = sus_pract_files, value = T)
dis_f_pract_files <- grep(pattern = "AF", x = dis_pract_files, value = T)
dis_m_pract_files <- grep(pattern = "AM", x = dis_pract_files, value = T)
afs_f_pract_files <- grep(pattern = "AM", x = afs_pract_files, value = T, invert = T)
afs_m_pract_files <- grep(pattern = "AM", x = afs_pract_files, value = T)

#read img in from final sess
sus_final %<-% lapply(list(sus_m_final_files, sus_f_final_files), function(x) readImage(x, 50, 68))
has_final %<-% lapply(list(has_m_final_files, has_f_final_files), function(x) readImage(x, 50, 68))
dis_final %<-% lapply(list(dis_m_final_files, dis_f_final_files), function(x) readImage(x, 50, 68))
ans_final %<-% lapply(list(ans_m_final_files, ans_f_final_files), function(x) readImage(x, 50, 68))

#read image in from practice sess males [[1]], females[[2]]
sus_pract %<-% lapply(list(sus_m_pract_files, sus_f_pract_files), function(x)readImage(x, 50,68))
ans_pract %<-% lapply(list(ans_m_pract_files, ans_f_pract_files), function(x)readImage(x, 50,68))
dis_pract %<-% lapply(list(dis_m_pract_files, dis_f_pract_files), function(x)readImage(x, 50,68))
afs_pract %<-% lapply(list(afs_m_pract_files, afs_f_pract_files), function(x)readImage(x, 50,68))

#transpose everything:
sus_final <- lapply(sus_final, transpose)
has_final <- lapply(has_final, transpose)
dis_final <- lapply(dis_final, transpose)
ans_final <- lapply(ans_final, transpose)
sus_pract <- lapply(sus_pract, transpose)
dis_pract <- lapply(dis_pract, transpose)
afs_pract <- lapply(afs_pract, transpose)

#1. return male disgust
#target: male dis:sus
#background: m-sus f-sus
tg1 <- rbind(dis_final[[1]], sus_final[[1]])
tg1 <- scale(tg1)
bg1 <- lapply(sus_pract, scale)
bg1_stack <- scale(do.call(rbind, sus_pract))

uca1 <- uca(tg1, bg1, nv = 5)
uca1_stack <- uca(tg1, bg1_stack, nv = 5)

#2. return male disgust
#target: male dis:sus
#background: m-sus, f-sus m-afs f-afs
bg2 <- lapply(c(sus_pract, afs_pract), scale)
bg2_stack <- scale(do.call(rbind,c(sus_pract, afs_pract)))

uca2 <- uca(tg1, bg2, nv = 5)
uca2_stack <- uca(tg1, bg2_stack, nv = 5)

#3. return female disgust
#target: female dis:sus
#background: f-sus, m-sus

tg3 <- rbind(dis_final[[2]], sus_final[[2]])
tg3 <- scale(tg3)
#bg3 = bg1
#bg3_stack = bg1_stack

uca3 <- uca(tg3, bg1, nv = 5)
uca3_stack <- uca(tg3, bg1_stack, nv = 5)

#4. return female disgust
#target: female dis:sus
#background: m-sus, f-sus m-afs f-afs

# tg4 = tg3
# bg4 = bg2

uca4 <- uca(tg3, bg2, nv = 5)
uca4_stack <- uca(tg3, bg2_stack, nv = 5)

#5. return male happy
#target: male has:sus
#background: m-sus f-sus
tg5 <- rbind(has_final[[1]], sus_final[[1]])
tg5 <- scale(tg5)
# bg5 = bg1

uca5 <- uca(tg5, bg1, nv = 5)
uca5_stack <- uca(tg5, bg1_stack, nv = 5)

#6. return male happy
#target: male has:sus
#background: m-sus, f-sus m-afs f-afs
#tg6 = tg5
#bg6  =  bg2

uca6 <- uca(tg5, bg2, nv = 5)
uca6_stack <- uca(tg5, bg2_stack, nv = 5)

#7. return: female happy
#target: female has:sus
#background: f-sus, m-sus
tg7 <- rbind(has_final[[2]], sus_final[[2]])
tg7 <- scale(tg7)
#bg7  =  bg1

uca7 <- uca(tg7, bg1, nv = 5)
uca7_stack <- uca(tg7, bg1_stack, nv = 5)

#8. return: female happy
#target: female has:sus
#background: m-sus, f-sus m-afs f-afs

#tg8 = tg7
#bg8 = bg2

uca8 <- uca(tg7, bg2, nv = 5)
uca8_stack <- uca(tg7, bg2_stack, nv = 5)

#9 return: male angry
#target: male ans:sus
#background: m-sus, m-afs

tg9 = rbind(sus_final[[1]], ans_final[[1]])
bg9 = lapply(list(sus_pract[[1]], afs_pract[[1]]), scale)
bg9_stack = scale(rbind(sus_pract[[1]], afs_pract[[1]]))
uca9 <- uca(tg9, bg9, nv = 5)
uca9_stack <- uca(tg9, bg9_stack, nv = 5)

#10 return: male angry
#target: male ans:sus
#background: m-sus, f-sus, m-afs, f-afs

#tg10=tg9
#bg10 = bg2
uca10 <- uca(tg9, bg2, nv = 5)
uca10_stack <- uca(tg9, bg2_stack, nv = 5)

# PCA of targets
#tg1, tg3, tg5, tg7,tg9
pca_1 <- eigs_sym(cov(tg1), k = 5, which = "LA")$vectors
pca_3 <- eigs_sym(cov(tg3), k = 5, which = "LA")$vectors
pca_5 <- eigs_sym(cov(tg5), k = 5, which = "LA")$vectors
pca_7 <- eigs_sym(cov(tg7), k = 5, which = "LA")$vectors
pca_9 <- eigs_sym(cov(tg9), k = 5, which = "LA")$vectors

# target 1 eigenfaces
tg1_list <- lapply(list(uca1, uca1_stack, uca2,uca2_stack), "[[", "vectors")
tg1_list <- lapply(tg1_list, function(mat) pos_corr(mat, pca_1))

plot_t1 <- rbind(toEigenfaces(pca_1, 50, 68, "PCA"),
                 toEigenfaces(tg1_list[[1]], 50, 68, "bg1"),
                 toEigenfaces(tg1_list[[2]], 50, 68, "bg1_stack"),
                 toEigenfaces(tg1_list[[3]], 50, 68, "bg2"),
                 toEigenfaces(tg1_list[[4]], 50, 68, "bg2_stack"))

plot_t1[, alpha := factor(alpha, levels = c("PCA", "bg1", "bg1_stack", "bg2","bg2_stack"))]
final_plot_t1 <- plotEigenfaces2(plot_t1, "")

ggsave("example/Faces/m_disgust_surprise.png",  width = 8, height = 11, units = "in")

# target 3 eigenfaces
tg3_list <- lapply(list(uca3, uca3_stack, uca4,uca4_stack), "[[", "vectors")
tg3_list <- lapply(tg3_list, function(mat) pos_corr(mat, pca_3))

plot_t3 <- rbind(toEigenfaces(pca_3, 50, 68, "PCA"),
                 toEigenfaces(tg3_list[[1]], 50, 68, "bg1"),
                 toEigenfaces(tg3_list[[2]], 50, 68, "bg1_stack"),
                 toEigenfaces(tg3_list[[3]], 50, 68, "bg2"),
                 toEigenfaces(tg3_list[[4]], 50, 68, "bg2_stack"))

plot_t3[, alpha := factor(alpha, levels = c("PCA", "bg1", "bg1_stack", "bg2","bg2_stack"))]
final_plot_t3 <- plotEigenfaces2(plot_t3, "")

ggsave("example/Faces/f_disgust_surprise.png",  width = 8, height = 11, units = "in")

# target 5 eigenfaces
tg5_list <- lapply(list(uca5, uca5_stack, uca6,uca6_stack), "[[", "vectors")
tg5_list <- lapply(tg5_list, function(mat) pos_corr(mat, pca_5))

plot_t5 <- rbind(toEigenfaces(pca_5, 50, 68, "PCA"),
                 toEigenfaces(tg5_list[[1]], 50, 68, "bg1"),
                 toEigenfaces(tg5_list[[2]], 50, 68, "bg1_stack"),
                 toEigenfaces(tg5_list[[3]], 50, 68, "bg2"),
                 toEigenfaces(tg5_list[[4]], 50, 68, "bg2_stack"))

plot_t5[, alpha := factor(alpha, levels = c("PCA", "bg1", "bg1_stack", "bg2","bg2_stack"))]
final_plot_t5 <- plotEigenfaces2(plot_t5, "")

ggsave("example/Faces/m_happy_surprise.png",  width = 8, height = 11, units = "in")

# target 7 eigenfaces
tg7_list <- lapply(list(uca7, uca7_stack, uca8,uca8_stack), "[[", "vectors")
tg7_list <- lapply(tg7_list, function(mat) pos_corr(mat, pca_7))

plot_t7 <- rbind(toEigenfaces(pca_7, 50, 68, "PCA"),
                 toEigenfaces(tg7_list[[1]], 50, 68, "bg1"),
                 toEigenfaces(tg7_list[[2]], 50, 68, "bg1_stack"),
                 toEigenfaces(tg7_list[[3]], 50, 68, "bg2"),
                 toEigenfaces(tg7_list[[4]], 50, 68, "bg2_stack"))

plot_t7[, alpha := factor(alpha, levels = c("PCA", "bg1", "bg1_stack", "bg2","bg2_stack"))]
final_plot_t7 <- plotEigenfaces2(plot_t7, "")
ggsave("example/Faces/f_happy_surprise.png",  width = 8, height = 11, units = "in")

# target 9 eigenfaces
tg9_list <- lapply(list(uca9, uca9_stack, uca10,uca10_stack), "[[", "vectors")
tg9_list <- lapply(tg9_list, function(mat) pos_corr(mat, pca_9))

plot_t9 <- rbind(toEigenfaces(pca_9, 50, 68, "PCA"),
                 toEigenfaces(tg9_list[[1]], 50, 68, "bg9"),
                 toEigenfaces(tg9_list[[2]], 50, 68, "bg9_stack"),
                 toEigenfaces(tg9_list[[3]], 50, 68, "bg2"),
                 toEigenfaces(tg9_list[[4]], 50, 68, "bg2_stack"))
 
plot_t9[, alpha := factor(alpha, levels = c("PCA", "bg9", "bg9_stack", "bg2","bg2_stack"))]
final_plot_t9 <- plotEigenfaces2(plot_t9, "")
ggsave("example/Faces/angry_surprise.png",  width = 8, height = 11, units = "in")

