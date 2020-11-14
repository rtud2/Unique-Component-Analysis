rm(list = ls())
Sys.setenv(MKL_DEBUG_CPU_TYPE=5)#, MKL_VERBOSE=1) #intel MKL settings. latter to hack amd processors
#setwd("/mnt/d/Residual-Dimension-Reduction")
setwd("D:/Residual-Dimension-Reduction/")
#install.packages("uca_0.13.zip", repos = NULL, type="source")
#install.packages("uca_0.13.tar.gz",type="source", repos = NULL) 
libraries <- c("data.table", "MASS", "ggplot2", "splines","gridExtra", "uca", "imager", "future", "RSpectra", "microbenchmark", "Rfast")

lapply(libraries, library, character.only = T)
plan(multiprocess)

## Helper functions 
source("functions/helper_functions.R")
kdef <- "faces/KDEF_and_AKDEF/KDEF/"
kdef_A <- "faces/KDEF_and_AKDEF/A_kdef_files/"
sus_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"SUS.JPG")))
sas_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"SAS.JPG")))
nes_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"NES.JPG")))
dis_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"DIS.JPG")))
ans_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"ANS.JPG")))
has_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"HAS.JPG")))

practice_list <- list(sus_pract_files, has_pract_files,sas_pract_files,nes_pract_files,dis_pract_files,ans_pract_files)

practice_male_list <- lapply(practice_list, function(zz){grep(pattern = "AM", x = zz, value = T)})
practice_female_list <- lapply(practice_list, function(zz){grep(pattern = "AF", x = zz, value = T)})

male_emotion_list %<-% lapply(practice_male_list, function(x) readImage(x, 50, 68))
female_emotion_list %<-% lapply(practice_female_list, function(x) readImage(x, 50, 68))

male_emotion_list <- lapply(male_emotion_list, transpose)
female_emotion_list <- lapply(female_emotion_list, transpose)

#pairwise female

var_sim <- function(d1, d2){
  # similarity between the stack and split matrices that make up the stack
  stk <- scale(rbind(d1,d2))
  cov_stk <- cov(stk)
  cov_d1 <- cov(scale(d1))
  cov_d2 <- cov(scale(d2))
  
  v <- eigs_sym(cov_stk, k = 1, which = "LA")$vectors
  stack_var <- transpose(v) %*% cov_stk %*%v
  d1_var <- transpose(v) %*% cov_d1 %*%v
  d2_var <- transpose(v) %*% cov_d2 %*%v
  cbind(stack_var, d1_var, d2_var)
}

inner_prod_sim <- function(d1, d2){
  # similarity between the stack and split matrices that make up the stack
  stk <- scale(rbind(d1,d2))
  cov_stk <- cov(stk)
  cov_d1 <- cov(scale(d1))
  cov_d2 <- cov(scale(d2))
  
  v_stack <- eigs_sym(cov_stk, k = 1, which = "LA")$vectors
  v_1 <- eigs_sym(cov_d1, k = 1, which = "LA")$vectors
  v_2 <- eigs_sym(cov_d2, k = 1, which = "LA")$vectors
  
  stack_1 <- crossprod(v_stack, v_1)
  stack_2 <- crossprod(v_stack, v_2)
  v_1_2 <- crossprod(v_1, v_2)
  cbind(stack_1, stack_2, v_1_2)
}

f_sus_comparison %<-% lapply(female_emotion_list[2:6], function(z) var_sim(d1 = female_emotion_list[[1]], z))
f_has_comparison %<-% lapply(female_emotion_list[3:6], function(z) var_sim(d1 = female_emotion_list[[2]], z))
f_sas_comparison %<-% lapply(female_emotion_list[4:6], function(z) var_sim(d1 = female_emotion_list[[3]], z))
f_nes_comparison %<-% lapply(female_emotion_list[5:6], function(z) var_sim(d1 = female_emotion_list[[4]], z))
f_dis_comparison %<-% var_sim(d1 = female_emotion_list[[6]], female_emotion_list[[5]])

comparison_num <- do.call(rbind,
        c(f_sus_comparison,
        f_has_comparison,
        f_sas_comparison,
        f_nes_comparison,
        list(f_dis_comparison)))

comparison_labs <- cbind(c(rep("sus", 5), rep("has", 4), rep("sas", 3), rep("nes", 2), "ans"),
      c("has","sas","nes","ans","dis",
        "sas","nes","ans","dis",
        "nes","ans","dis",
        "ans","dis",
        "dis"))

background_var <- data.table(comparison_labs, comparison_num)
setnames(background_var, c("bg1","bg2","Stack_var","bg1_var", "bg2_var"))

background_var[Stack_var > 0.5*(bg1_var +  bg2_var)]
background_var[, diff := abs(bg1_var - bg2_var)] #if diff is large means stack is overrun by one group

background_var[order(diff)]
# bg1 bg2 Stack_var  bg1_var  bg2_var       diff
# 1: has ans  797.7054 821.4231 825.0023   3.579124
# 2: sas dis  763.4032 803.3573 772.8698  30.487463
# 3: has sas  783.9296 830.6244 796.9264  33.698050
# 4: sas ans  800.8645 796.8466 832.2659  35.419341
# 5: sus nes  872.2886 915.5720 875.5956  39.976461
# 6: nes ans  793.3335 855.7559 814.1901  41.565779
# 7: has nes  825.7674 831.2871 876.8953  45.608142
# 8: has dis  747.6688 819.9435 762.3158  57.627715
# 9: ans dis  766.7865 760.2942 824.7095  64.415299
# 10: sus ans  817.5896 890.6497 807.1156  83.534117
# 11: sas nes  811.0090 793.3120 878.0899  84.777978
# 12: sus has  841.1609 910.9744 826.0987  84.875680
# 13: nes dis  770.3364 863.9260 754.1740 109.751996
# 14: sus sas  829.5219 911.0035 780.1336 130.869922
# 15: sus dis  786.1118 896.8199 734.6812 162.138647

# sus & dis
# sus & sas
# nes & dis


# INNER PRODUCT
f_sus_inner %<-% lapply(female_emotion_list[2:6], function(z) inner_prod_sim(d1 = female_emotion_list[[1]], z))
f_has_inner %<-% lapply(female_emotion_list[3:6], function(z) inner_prod_sim(d1 = female_emotion_list[[2]], z))
f_sas_inner %<-% lapply(female_emotion_list[4:6], function(z) inner_prod_sim(d1 = female_emotion_list[[3]], z))
f_nes_inner %<-% lapply(female_emotion_list[5:6], function(z) inner_prod_sim(d1 = female_emotion_list[[4]], z))
f_dis_inner %<-% inner_prod_sim(d1 = female_emotion_list[[6]], female_emotion_list[[5]])

inner_num <- do.call(rbind,
        c(f_sus_inner,
        f_has_inner,
        f_sas_inner,
        f_nes_inner,
        list(f_dis_inner)))

inner_labs <- cbind(c(rep("sus", 5), rep("has", 4), rep("sas", 3), rep("nes", 2), "ans"),
      c("has","sas","nes","ans","dis",
        "sas","nes","ans","dis",
        "nes","ans","dis",
        "ans","dis",
        "dis"))

background_inner <- data.table(inner_labs, inner_num)
setnames(background_inner, c("bg1","bg2","Stack_1","Stack_2", "bg1_bg2"))
max_inner <- background_inner[, .(max = max(Stack_1, Stack_2, bg1_bg2)), by = c("bg1","bg2")]
max_inner[order(max)]
# bg1 bg2       max
# 1: sus ans 0.9691885
# 2: nes ans 0.9727309
# 3: nes dis 0.9791950
# 4: sus dis 0.9795358
# 5: ans dis 0.9812969
# 6: has ans 0.9814274
# 7: has dis 0.9838760
# 8: sas nes 0.9862520
# 9: sus sas 0.9864356
# 10: sas ans 0.9868335
# 11: has sas 0.9868691
# 12: sus has 0.9872419
# 13: has nes 0.9881727
# 14: sus nes 0.9905884
# 15: sas dis 0.9909845
# 
background_inner[order(bg1_bg2)]
# bg1 bg2   Stack_1   Stack_2   bg1_bg2
# 1: sus ans 0.9691885 0.9665842 0.8800309
# 2: nes ans 0.9659761 0.9727309 0.8864153
# 3: sus dis 0.9795358 0.9586367 0.8959284
# 4: nes dis 0.9791950 0.9742804 0.9166890
# 5: sus sas 0.9864356 0.9676066 0.9167276
# 6: has ans 0.9766456 0.9814274 0.9220293
# 7: sas nes 0.9777867 0.9862520 0.9324228
# 8: ans dis 0.9810512 0.9812969 0.9332532
# 9: has dis 0.9803653 0.9838760 0.9366546
# 10: sas ans 0.9816926 0.9868335 0.9409767
# 11: has sas 0.9868691 0.9824078 0.9441928
# 12: sus has 0.9872419 0.9843728 0.9531639
# 13: sus nes 0.9905884 0.9867179 0.9575075
# 14: sas dis 0.9869130 0.9909845 0.9597885
# 15: has nes 0.9881727 0.9867766 0.9606486
# read in final pictures
sus_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"SUS.JPG")))
has_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"HAS.JPG")))
sas_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"SAS.JPG")))
nes_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"NES.JPG")))
dis_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"DIS.JPG")))
ans_final_files <- paste0(kdef, sapply(list.files(kdef), function(zz) paste0(zz,"/",zz,"ANS.JPG")))

final_list <- list(sus_final_files, has_final_files,sas_final_files,nes_final_files,dis_final_files,ans_final_files)

# final_male_list <- lapply(final_list, function(zz){grep(pattern = "AM", x = zz, value = T)})
# final_male_emotion_list %<-% lapply(final_male_list, function(x) readImage(x, 50, 68))
# final_male_emotion_list <- lapply(final_male_emotion_list, transpose)
final_female_list <- lapply(final_list, function(zz){grep(pattern = "BF", x = zz, value = T)})
final_female_emotion_list %<-% lapply(final_female_list, function(x) readImage(x, 50, 68))
final_female_emotion_list <- lapply(final_female_emotion_list, transpose)


targ_var_sim <- function(targ, d1, d2){
  # similarity between the stack and split matrices that make up the stack
  # to the target data
  stk <- scale(rbind(d1,d2))
  cov_stk <- cov(stk)
  cov_d1 <- cov(scale(d1))
  cov_d2 <- cov(scale(d2))
  cov_targ <- cov(scale(targ))
  
  v <- eigs_sym(cov_stk, k = 1, which = "LA")$vectors
  targ_var <- transpose(v) %*% cov_targ %*% v
  stack_var <- transpose(v) %*% cov_stk %*%v
  d1_var <- transpose(v) %*% cov_d1 %*%v
  d2_var <- transpose(v) %*% cov_d2 %*%v
  cbind(targ_var, stack_var, d1_var, d2_var)
}
targ_inner_prod_sim <- function(targ, d1, d2){
  # similarity between the stack and split matrices that make up the stack
  stk <- scale(rbind(d1,d2))
  cov_stk <- cov(stk)
  cov_d1 <- cov(scale(d1))
  cov_d2 <- cov(scale(d2))
  cov_targ <- cov(scale(targ))
  
  v_targ <- eigs_sym(cov_targ, k = 1, which = "LA")$vectors
  v_stack <- eigs_sym(cov_stk, k = 1, which = "LA")$vectors
  v_1 <- eigs_sym(cov_d1, k = 1, which = "LA")$vectors
  v_2 <- eigs_sym(cov_d2, k = 1, which = "LA")$vectors
  
  targ_stack <- crossprod(v_targ, v_stack)
  targ_1 <- crossprod(v_targ, v_1)
  targ_2 <- crossprod(v_targ, v_2)
  cbind(targ_stack, targ_1, targ_2)
}

#sus ans
#order: has,sas,nes,dis,ans
sus_ans <- do.call(rbind, lapply(final_female_emotion_list[c(-1,-6)],
       function(zz){ targ_inner_prod_sim(zz, female_emotion_list[[1]], female_emotion_list[[6]])}))
sus_ans_inner <- data.table(c("HAS","SAS","NES","DIS") , sus_ans)
setnames(sus_ans_inner, c("Target","Inner_Stack","Inner_1","Inner_2"))
# Target Inner_Stack   Inner_1   Inner_2
# 1:    HAS   0.9565713 0.9435381 0.9333191
# 2:    SAS   0.9538911 0.9236241 0.9466494
# 3:    NES   0.9604394 0.9468185 0.9227120
# 4:    DIS   0.9518510 0.9151435 0.9597355

# nes ans 
#order: sus,has,sas,dis
nes_ans <- do.call(rbind, lapply(final_female_emotion_list[c(-4,-6)],
       function(zz){ targ_inner_prod_sim(zz, female_emotion_list[[4]], female_emotion_list[[6]])}))
nes_ans_inner <- data.table(c("SUS","HAS","SAS","DIS") , sus_ans)
setnames(nes_ans_inner, c("Target","Inner_Stack","Inner_1","Inner_2"))
# Target Inner_Stack   Inner_1   Inner_2
# 1:    SUS   0.9565713 0.9435381 0.9333191
# 2:    HAS   0.9538911 0.9236241 0.9466494
# 3:    SAS   0.9604394 0.9468185 0.9227120
# 4:    DIS   0.9518510 0.9151435 0.9597355

# nes dis 
#order: sus,has,sas,ans
nes_dis <- do.call(rbind, lapply(final_female_emotion_list[c(-4,-5)],
       function(zz){ targ_inner_prod_sim(zz, female_emotion_list[[4]], female_emotion_list[[5]])}))
nes_dis_inner <- data.table(c("SUS","HAS","SAS","ANS") , sus_ans)
setnames(nes_dis_inner, c("Target","Inner_Stack","Inner_1","Inner_2"))
# Target Inner_Stack   Inner_1   Inner_2
# 1:    SUS   0.9565713 0.9435381 0.9333191
# 2:    HAS   0.9538911 0.9236241 0.9466494
# 3:    SAS   0.9604394 0.9468185 0.9227120
# 4:    ANS   0.9518510 0.9151435 0.9597355

#sus dis
f_sus_dis <- do.call(rbind, final_female_emotion_list[c(1,5)])
f_not_sus_dis <- final_female_emotion_list[-c(1,5)]

f_not_sus_dis_cat <- lapply(f_not_sus_dis, function(zz){rbind(f_sus_dis, zz)})

sus_dis %<-% lapply(f_not_sus_dis_cat,
                  function(z) targ_var_sim(targ = z,
                                           d1 = female_emotion_list[[1]],
                                           d2 = female_emotion_list[[5]]))
#order idx: has,sas,nes,ans

#sus sas
f_sus_sas <- do.call(rbind, final_female_emotion_list[c(1,3)])
f_not_sus_sas <- final_female_emotion_list[-c(1,3)]

f_not_sus_sas_cat <- lapply(f_not_sus_sas, function(zz){rbind(f_sus_sas, zz)})

sus_sas %<-% lapply(f_not_sus_sas_cat,
                  function(z) targ_var_sim(targ = z,
                                           d1 = female_emotion_list[[1]],
                                           d2 = female_emotion_list[[3]]))
#order idx: has,nes,dis,ans

#nes dis
f_nes_dis <- do.call(rbind, final_female_emotion_list[c(4,5)])
f_not_nes_dis <- final_female_emotion_list[-c(4,5)]

f_not_nes_dis_cat <- lapply(f_not_nes_dis, function(zz){rbind(f_nes_dis, zz)})

nes_dis %<-% lapply(f_not_nes_dis_cat,
                  function(z) targ_var_sim(targ = z,
                                           d1 = female_emotion_list[[4]],
                                           d2 = female_emotion_list[[5]]))

#order idx: sus,has,sas,ans

comparison_target <- rbind(
data.table(c("has","sas","nes","ans"),"sus_dis", do.call(rbind, sus_dis)),
data.table(c("has","sas","nes","ans"),"sus_sas", do.call(rbind, sus_sas)),
data.table(c("sus","has","sas","ans"),"nes_dis", do.call(rbind, nes_dis)))

setnames(comparison_target, c("Target","Background", "Target.Variation", "Stack.Variation","Bg1.Variation","Bg2.Variation"))

comparison_target[,.SD, keyby = c("Background","Target.Variation")]
# Background Target.Variation Target Stack.Variation Bg1.Variation Bg2.Variation
# 1:    nes_dis         702.5117    ans        793.3335      855.7559      814.1901
# 2:    nes_dis         745.8252    sas        793.3335      855.7559      814.1901
# 3:    nes_dis         761.1969    has        793.3335      855.7559      814.1901
# 4:    nes_dis         767.0170    sus        793.3335      855.7559      814.1901

# 5:    sus_dis         715.5389    ans        817.5896      890.6497      807.1156
# 6:    sus_dis         762.3043    sas        817.5896      890.6497      807.1156
# 7:    sus_dis         763.4855    nes        817.5896      890.6497      807.1156
# 8:    sus_dis         776.6644    has        817.5896      890.6497      807.1156

# 9:    sus_sas         719.0278    ans        829.5219      911.0035      780.1336
# 10:    sus_sas         763.0583    nes        829.5219      911.0035      780.1336
# 11:    sus_sas         767.5979    has        829.5219      911.0035      780.1336
# 12:    sus_sas         773.3774   sas        829.5219      911.0035      780.1336
#

# combos that stack does the worst in target, but also washes out one of the signals
# final_list <- list(sus_final_files, has_final_files,sas_final_files,nes_final_files,dis_final_files,ans_final_files)
# practice_list <- list(sus_pract_files, has_pract_files,sas_pract_files,nes_pract_files,dis_pract_files,ans_pract_files)

#nes-dis-ans 
#target: final_female_emotion_list[6] (ANS)
#background: female_emotion_list[c(4,5)] (NES,DIS)
nes_bg <- scale(female_emotion_list[[4]])
dis_bg <- scale(female_emotion_list[[5]])
nes_dis_stack <- scale(do.call(rbind, female_emotion_list[c(4,5)]))
nes_dis_ans <- scale(do.call(rbind, final_female_emotion_list[c(4,5,6)]))

uca1_nes_ans <- uca(nes_dis_ans, nes_bg, nv = 5)
uca1_dis_ans <- uca(nes_dis_ans, dis_bg, nv = 5)
uca1_nes_dis_ans <- uca(nes_dis_ans, list(nes_bg, dis_bg), nv = 5)
uca1_nes_dis_ans_stk <- uca(nes_dis_ans, nes_dis_stack, nv = 5)

nes_dis_stack_cov <- cov(nes_dis_stack)
diag(nes_dis_stack_cov) <- diag(nes_dis_stack_cov) + 0.0001
cpcapp1 <- eigs_sym(solve(nes_dis_stack_cov) %*% cov(nes_dis_ans), k = 5 , which = "LA")$vectors

pca_nes_dis_ans <- eigs_sym(cov(nes_dis_ans), k = 5, which = "LA")$vectors
pca_ans <- eigs_sym(cov(final_female_emotion_list[[6]]), k = 5, which = "LA")$vectors

nes_dis_ans_cov_list <- lapply(final_female_emotion_list[c(4,5,6)], function(z)cov(scale(z)))

#sus-dis-ans
#target: final_female_emotion_list[6] (ANS)
#background: female_emotion_list[c(1,5)] (SUS,DIS)
sus_bg <- scale(female_emotion_list[[1]])
sus_dis_stack <- scale(do.call(rbind, female_emotion_list[c(1,5)]))
sus_dis_ans <- scale(do.call(rbind, final_female_emotion_list[c(1,5,6)]))

uca2_sus_ans <- uca(sus_dis_ans, sus_bg, nv = 5)
uca2_dis_ans <- uca(sus_dis_ans, dis_bg, nv = 5)
uca2_sus_dis_ans <- uca(sus_dis_ans, list(sus_bg, dis_bg), nv = 5)
uca2_sus_dis_ans_stk <- uca(sus_dis_ans, sus_dis_stack, nv = 5)

sus_dis_stack_cov <- cov(sus_dis_stack)
diag(sus_dis_stack_cov) <- diag(sus_dis_stack_cov) + 0.0001
cpcapp2 <- eigs_sym(solve(sus_dis_stack_cov) %*% cov(sus_dis_ans), k = 5 , which = "LA")$vectors

pca_sus_dis_ans <- eigs_sym(cov(sus_dis_ans), k = 5, which = "LA")$vectors

sus_dis_ans_cov_list <- lapply(final_female_emotion_list[c(1,5,6)], function(z)cov(scale(z)))

#sus-sas-ans
#target: final_female_emotion_list[6] (ANS)
#background: female_emotion_list[c(1,3)] (SUS,SAS)
sas_bg <- scale(female_emotion_list[[3]])
sus_sas_stack <- scale(do.call(rbind, female_emotion_list[c(1,3)]))
sus_sas_ans <- scale(do.call(rbind, final_female_emotion_list[c(1,3,6)]))

uca3_sus_ans <- uca(sus_sas_ans, sus_bg, nv = 5)
uca3_sas_ans <- uca(sus_sas_ans, sas_bg, nv = 5)
uca3_sus_sas_ans <- uca(sus_sas_ans, list(sus_bg, sas_bg), nv = 5)
uca3_sus_sas_ans_stk <- uca(sus_sas_ans, sus_sas_stack, nv = 5)

sus_sas_stack_cov <- cov(sus_sas_stack)
diag(sus_sas_stack_cov) <- diag(sus_sas_stack_cov) + 0.0001
cpcapp3 <- eigs_sym(solve(sus_sas_stack_cov) %*% cov(sus_sas_ans), k = 5 , which = "LA")$vectors

pca_sus_sas_ans <- eigs_sym(cov(sus_sas_ans), k = 5, which = "LA")$vectors
sus_sas_ans_cov_list <- lapply(final_female_emotion_list[c(1,3,6)], function(z)cov(scale(z)))

#ans-dis-sus
#target: final_female_emotion_list[1] (SUS)
#background: female_emotion_list[c(5,6)] (DIS,ANS)
ans_bg <- scale(female_emotion_list[[6]])
ans_dis_stack <- scale(do.call(rbind, female_emotion_list[c(5,6)]))

uca4_ans_sus <- uca(sus_dis_ans, ans_bg, nv = 5)
uca4_dis_sus <- uca(sus_dis_ans, dis_bg, nv = 5)
uca4_ans_dis_sus <- uca(sus_dis_ans, list(ans_bg, dis_bg), nv = 5)
uca4_ans_dis_sus_stk <- uca(sus_dis_ans, ans_dis_stack, nv = 5)

ans_dis_stack_cov <- cov(ans_dis_stack)
diag(ans_dis_stack_cov) <- diag(ans_dis_stack_cov) + 0.0001
cpcapp4 <- eigs_sym(solve(ans_dis_stack_cov) %*% cov(sus_dis_ans), k = 5 , which = "LA")$vectors

pca_sus <- eigs_sym(cov(final_female_emotion_list[[1]]), k = 5, which = "LA")$vectors


getFaceProjection=function(idx1, idx2){
  emo_names <- c("SUS","HAS", "SAS", "NES", "DIS","ANS")
  # function to cat class names onto output...
  bg_idx <- c(idx1,idx2)
  class_names = lapply(emo_names[-bg_idx], function(others){rep(c(others,emo_names[bg_idx]), each = 35)})  
  
  uca_bg1_bg2 <- lapply(final_female_emotion_list[-c(idx1,idx2)], function(targ){
    tmp_targ <- scale(do.call(rbind, c(list(targ), final_female_emotion_list[c(idx1,idx2)] )))
    tmp_bg1 <- scale(final_female_emotion_list[[idx1]])
    tmp_bg2 <- scale(final_female_emotion_list[[idx2]])
    tmp_stack <- scale(do.call(rbind, final_female_emotion_list[c(idx1,idx2)]))
    cov_tmp_stack <- cov(tmp_stack)
    diag(cov_tmp_stack) <- diag(cov_tmp_stack) + 1e-3
    
    targ_pca <- eigs_sym(cov(scale(targ)), k = 5, which = "LA")$vectors
    
    tmp_pca <- eigs_sym(cov(tmp_targ), k = 5, which = "LA")$vectors
    tmp_cpcapp <- eigs_sym(solve(cov_tmp_stack) %*% cov(tmp_targ), k = 5 , which = "LA")$vectors
    uca_bg1 <- uca(tmp_targ, tmp_bg1, nv = 5)$vectors
    uca_bg2 <- uca(tmp_targ, tmp_bg2, nv = 5)$vectors
    uca_bg12 <- uca(tmp_targ, list(tmp_bg1, tmp_bg2), nv = 5)$vectors
    uca_bg12_stk <- uca(tmp_targ, tmp_stack, nv = 5)$vectors
    out_method <- c("PCA","cPCA++",emo_names[bg_idx], "Split", "Stack")
    
    out <- list(tmp_pca, tmp_cpcapp, uca_bg1, uca_bg2, uca_bg12, uca_bg12_stk )
    out <- lapply(out, function(mat) pos_corr(mat, targ_pca))
    
    projections <- do.call(rbind, lapply(seq_along(out), function(vec) data.table(tmp_targ %*% out[[vec]],
                                                                                  "method" = out_method[vec]  )))
    out <- data.table(do.call(rbind, out), "method" = factor(rep(out_method, each = 3400),out_method))
    list(out, projections) 
  })
  uca_plot <- lapply(uca_bg1_bg2 , "[[", 1)
  uca_proj <- lapply(uca_bg1_bg2, "[[", 2) #stores the projected eigenfaces on to 2-dims
  uca_proj <- Map(cbind, uca_proj, class_names)
  uca_proj <- lapply(uca_proj, function(dat) setnames(dat, c(paste0("PC",1:5), "Method","Class")))
  list("Plot" = uca_plot, "Proj" = uca_proj) 
}

#final_female_emotion_list: 1. sus, 2. has, 3. sas, 4. nes, 5. dis, 6. ans
#sus ans
#list order: has,sas,nes,dis,ans
#within list order: PCA, cPCA++, SUS, ANS, Split, Stack
uca_sus_ans <- getFaceProjection(1,6)
# nes ans  (4,6)
#order: sus,has,sas,dis
#classes <- c("SUS","HAS","SAS","DIS")
uca_nes_ans <- getFaceProjection(4,6)
# nes dis  (4,5)
#order: sus,has,sas,ans
uca_nes_dis <- getFaceProjection(4,5)


#plot the most contrastive...
uca1_list <- lapply(list(uca1_nes_ans, uca1_dis_ans, uca1_nes_dis_ans,uca1_nes_dis_ans_stk), "[[", "vectors")
uca1_list <- lapply(uca1_list, function(mat) pos_corr(mat, pca_ans))
pca_nes_dis_ans <- pos_corr(pca_nes_dis_ans,pca_ans)

plot_uca1 <- rbind(toEigenfaces(pca_nes_dis_ans, 50, 68, "PCA"),
                 toEigenfaces(uca1_list[[1]], 50, 68, "Neutral"),
                 toEigenfaces(uca1_list[[2]], 50, 68, "Disgust"),
                 toEigenfaces(uca1_list[[3]], 50, 68, "Split"),
                 toEigenfaces(uca1_list[[4]], 50, 68, "Stack"))

plot_uca1[, alpha := factor(alpha, levels = c("PCA", "Neutral", "Disgust", "Split","Stack"))]
final_plot_uca1 <- plotEigenfaces2(plot_uca1, "")

ggsave("example/Faces/F_Neutral_Disgust_Angry.png",  width = 8, height = 11, units = "in")

uca2_list <- lapply(list(uca2_sus_ans, uca2_dis_ans, uca2_sus_dis_ans,uca2_sus_dis_ans_stk), "[[", "vectors")
uca2_list <- lapply(uca2_list, function(mat) pos_corr(mat, pca_ans))
pca_sus_dis_ans <- pos_corr(pca_sus_dis_ans, pca_ans)
plot_uca2 <- rbind(toEigenfaces(pca_sus_dis_ans, 50, 68, "PCA"),
                 toEigenfaces(uca2_list[[1]], 50, 68, "Surprise"),
                 toEigenfaces(uca2_list[[2]], 50, 68, "Disgust"),
                 toEigenfaces(uca2_list[[3]], 50, 68, "Split"),
                 toEigenfaces(uca2_list[[4]], 50, 68, "Stack"))

plot_uca2[, alpha := factor(alpha, levels = c("PCA", "Surprise", "Disgust", "Split","Stack"))]
final_plot_uca2 <- plotEigenfaces2(plot_uca2, "")

ggsave("example/Faces/F_Surprise_Disgust_Angry.png",  width = 8, height = 11, units = "in")


uca3_list <- lapply(list(uca3_sus_ans, uca3_sas_ans, uca3_sus_sas_ans,uca3_sus_sas_ans_stk), "[[", "vectors")
uca3_list <- lapply(uca3_list, function(mat) pos_corr(mat, pca_ans))
pca_sus_sas_ans <- pos_corr(pca_sus_sas_ans, pca_ans)
plot_uca3 <- rbind(toEigenfaces(pca_sus_sas_ans, 50, 68, "PCA"),
                 toEigenfaces(uca3_list[[1]], 50, 68, "Surprise"),
                 toEigenfaces(uca3_list[[2]], 50, 68, "Sad"),
                 toEigenfaces(uca3_list[[3]], 50, 68, "Split"),
                 toEigenfaces(uca3_list[[4]], 50, 68, "Stack"))

plot_uca3[, alpha := factor(alpha, levels = c("PCA", "Surprise", "Sad", "Split","Stack"))]
final_plot_uca3 <- plotEigenfaces2(plot_uca3, "")

ggsave("example/Faces/F_Surprise_Sad_Angry.png",  width = 8, height = 11, units = "in")


uca4_list <- lapply(list(uca4_ans_sus, uca4_dis_sus, uca4_ans_dis_sus,uca4_ans_dis_sus_stk), "[[", "vectors")
uca4_list <- lapply(uca4_list, function(mat) pos_corr(mat, pca_ans))
pca_sus_dis_ans <- pos_corr(pca_sus_dis_ans, pca_sus)
plot_uca4 <- rbind(toEigenfaces(pca_sus_dis_ans, 50, 68, "PCA"),
                 toEigenfaces(uca4_list[[1]], 50, 68, "Angry"),
                 toEigenfaces(uca4_list[[2]], 50, 68, "Disgust"),
                 toEigenfaces(uca4_list[[3]], 50, 68, "Split"),
                 toEigenfaces(uca4_list[[4]], 50, 68, "Stack"))

plot_uca4[, alpha := factor(alpha, levels = c("PCA", "Angry", "Disgust", "Split","Stack"))]
final_plot_uca4 <- plotEigenfaces2(plot_uca4, "")

ggsave("example/Faces/F_Angry_Disgust_Surprise.png",  width = 8, height = 11, units = "in")


#sus ans
#list order: has,sas,nes,dis
#within list order: PCA, cPCA++, SUS, ANS, Split, Stack
sus_ans_emo_list <- c("Happy","Sad","Neutral","Disgust")
uca_sus_ans_plot <- uca_sus_ans[[1]]

sus_ans_plot <- lapply(seq_along(uca_sus_ans_plot), function(idx){
  do.call(rbind, lapply(uca_sus_ans_plot[[idx]][,unique(method)], function(inner){
    toEigenfaces(data.matrix(uca_sus_ans_plot[[idx]][method == inner, .SD, .SDcols = -"method"]), 50, 68, inner)}
    ))
  } )
#plot eigenfaces
lapply(seq_along(sus_ans_emo_list), function(zz){
  plotEigenfaces2(sus_ans_plot[[zz]],"")   
  ggsave(paste0("example/Faces/Surprise_Angry/F_Surprise_Angry_",sus_ans_emo_list[zz],".png"),  width = 8, height = 11, units = "in")
  })
#plot projected faces
lapply(seq_along(uca_sus_ans[[2]]), function(idx){
  tmp_dat <- uca_sus_ans[[2]][[idx]]
  tmp_dat[,Method := factor(Method, levels = c("PCA","cPCA++","SUS","ANS","Split","Stack"))]
  ggplot(data = tmp_dat, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(title = paste0("Female Surprise, Angry, ", sus_ans_emo_list[idx]), x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave(paste0("example/Faces/Surprise_Angry/F_Surprise_Angry_",sus_ans_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
    }        
  )

# nes ans  (4,6)
#order: sus,has,sas,dis
nes_ans_emo_list <- c("Surprise","Happy","Sad","Disgust")
uca_nes_ans_plot <- uca_nes_ans[[1]]

nes_ans_plot <- lapply(seq_along(uca_nes_ans_plot), function(idx){
  do.call(rbind, lapply(uca_nes_ans_plot[[idx]][,unique(method)], function(inner){
    toEigenfaces(data.matrix(uca_nes_ans_plot[[idx]][method == inner, .SD, .SDcols = -"method"]), 50, 68, inner)}
    ))
  } )
#plot eigenfaces
lapply(seq_along(nes_ans_emo_list), function(zz){
  plotEigenfaces2(nes_ans_plot[[zz]],"")   
  ggsave(paste0("example/Faces/Neutral_Angry/F_Neutral_Angry_",nes_ans_emo_list[zz],".png"),  width = 8, height = 11, units = "in")
  })
#plot projected faces
lapply(seq_along(uca_nes_ans[[2]]), function(idx){
  tmp_dat <- uca_nes_ans[[2]][[idx]]
  tmp_dat[,Method := factor(Method, levels = c("PCA","cPCA++","NES","ANS","Split","Stack"))]
  ggplot(data = tmp_dat, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(title = paste0("Female Neutral, Angry, ", nes_ans_emo_list[idx]), x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave(paste0("example/Faces/Neutral_Angry/F_Surprise_Angry_",nes_ans_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
    }        
  )

# nes dis  (4,5)
#order: sus,has,sas,ans
nes_dis_emo_list <- c("Surprise","Happy","Sad","Angry")
uca_nes_dis_plot <- uca_nes_dis[[1]]

nes_dis_plot <- lapply(seq_along(uca_nes_dis_plot), function(idx){
  do.call(rbind, lapply(uca_nes_dis_plot[[idx]][,unique(method)], function(inner){
    toEigenfaces(data.matrix(uca_nes_dis_plot[[idx]][method == inner, .SD, .SDcols = -"method"]), 50, 68, inner)}
    ))
  } )
#plot eigenfaces
lapply(seq_along(nes_dis_emo_list), function(zz){
  plotEigenfaces2(nes_dis_plot[[zz]],"")   
  ggsave(paste0("example/Faces/Neutral_Disgust/F_Neutral_Disgust_",nes_dis_emo_list[zz],".png"),  width = 8, height = 11, units = "in")
  })
#plot projected faces
lapply(seq_along(uca_nes_dis[[2]]), function(idx){
  tmp_dat <- uca_nes_dis[[2]][[idx]]
  tmp_dat[,Method := factor(Method, levels = c("PCA","cPCA++","NES","DIS","Split","Stack"))]
  ggplot(data = tmp_dat, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(title = paste0("Female Neutral, Disgust, ", nes_dis_emo_list[idx]), x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave(paste0("example/Faces/Neutral_Disgust/F_Surprise_Disgust_",nes_dis_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
    }        
  )

#final session variance explained (these should not be controlled)
nes_dis_ans_var_exp <- do.call(rbind,
        lapply(uca1_list, function(V){
  sapply(nes_dis_ans_cov_list,
         function(C) crossprod(V[,1,drop = F], C %*% V[,1,drop = F]))
    }))
colnames(nes_dis_ans_var_exp) <- c("Neutral","Disgust","Angry")
nes_dis_ans_var_exp <- data.table(nes_dis_ans_var_exp)
nes_dis_ans_var_exp[, Background :=  c("NES","DIS","Split","Stack")]


sus_dis_ans_var_exp <- do.call(rbind,
        lapply(uca2_list, function(V){
  sapply(sus_dis_ans_cov_list,
         function(C) crossprod(V[,1,drop = F], C %*% V[,1,drop = F]))
    }))

colnames(sus_dis_ans_var_exp) <- c("Surprise","Disgust","Angry")
sus_dis_ans_var_exp <- data.table(sus_dis_ans_var_exp)
sus_dis_ans_var_exp[, Background :=  c("SUS","DIS","Split","Stack")]

sus_sas_ans_var_exp <- do.call(rbind,
        lapply(uca3_list, function(V){
  sapply(sus_sas_ans_cov_list,
         function(C) crossprod(V[,1,drop = F], C %*% V[,1,drop = F]))
    }))

colnames(sus_sas_ans_var_exp) <- c("Surprise","Sad","Angry")
sus_sas_ans_var_exp <- data.table(sus_sas_ans_var_exp)
sus_sas_ans_var_exp[, Background :=  c("SUS","SAS","Split","Stack")]

ans_dis_sus_var_exp <- do.call(rbind,
        lapply(uca4_list, function(V){
  sapply(sus_dis_ans_cov_list,
         function(C) crossprod(V[,1,drop = F], C %*% V[,1,drop = F]))
    }))

colnames(ans_dis_sus_var_exp) <- c("Surprise","Disgust","Angry")
ans_dis_sus_var_exp <- data.table(sus_dis_ans_var_exp)
ans_dis_sus_var_exp[, Background :=  c("ANS","DIS","Split","Stack")]

fwrite(nes_dis_ans_var_exp, file = "example/Faces/NES_DIS_ANS_Var.csv")
fwrite(sus_dis_ans_var_exp, file = "example/Faces/SUS_DIS_ANS_Var.csv")
fwrite(sus_sas_ans_var_exp, file = "example/Faces/SUS_SAS_ANS_Var.csv")
fwrite(ans_dis_sus_var_exp, file = "example/Faces/ANS_DIS_SUS_Var.csv")

#classification results on 2 dimensions:
uca1_vectors <- c(list(pca_nes_dis_ans),list(cpcapp1), uca1_list)
uca1_methods <- c("PCA","cPCA++","Neutral","Disgust","Split","Stack")

uca1_plot <- do.call(rbind,
lapply(seq_along(uca1_vectors), function(V) {
  res <- data.table(nes_dis_ans %*% uca1_vectors[[V]][,1:2], rep(c("NES","DIS","ANS"), rep(35,3)),uca1_methods[V])
  setnames(res, c("PC1", "PC2", "Class","Method"))
  }))
uca1_plot[, Method := factor(Method, levels = c("PCA","cPCA++","Neutral","Disgust","Split","Stack"))]
ggplot(data = uca1_plot, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(title = "Female, Neutral, Disgust, Angry", x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("example/Faces/F_Neutral_Disgust_Angry_projected.png",  width = 11, height = 8, units = "in")

uca2_vectors <- c(list(pca_sus_dis_ans),list(cpcapp2), uca2_list)
uca2_methods <- c("PCA","cPCA++","Surprise","Disgust","Split","Stack")

uca2_plot <- do.call(rbind,
lapply(seq_along(uca2_vectors), function(V) {
  res <- data.table(sus_dis_ans %*% uca2_vectors[[V]][,1:2], rep(c("SUS","DIS","ANS"), rep(35,3)),uca2_methods[V])
  setnames(res, c("PC1", "PC2", "Class","Method"))
  }))
uca2_plot[,Method := factor(Method, levels = c("PCA","cPCA++","Surprise", "Disgust","Split","Stack"))]
ggplot(data = uca2_plot, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(title = "Female Surprise, Disgust, Angry", x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("example/Faces/F_Surprise_Disgust_Angry_projected.png",  width = 11, height = 8, units = "in")


uca3_vectors <- c(list(pca_sus_sas_ans),list(cpcapp3), uca3_list)
uca3_methods <- c("PCA","cPCA++","Surprise","Sad","Split","Stack")

uca3_plot <- do.call(rbind,
lapply(seq_along(uca3_vectors), function(V) {
  res <- data.table(sus_sas_ans %*% uca3_vectors[[V]][,1:2], rep(c("SUS","SAS","ANS"), rep(35,3)),uca3_methods[V])
  setnames(res, c("PC1", "PC2", "Class","Method"))
  }))
uca3_plot[,Method := factor(Method, levels = c("PCA","cPCA++","Surprise", "Sad","Split","Stack"))]
ggplot(data = uca3_plot, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(title = "Female Surprise, Sad, Angry", x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("example/Faces/F_Surprise_Sad_Angry_projected.png",  width = 11, height = 8, units = "in")
 
uca4_vectors <- c(list(pca_sus_dis_ans),list(cpcapp4), uca4_list)
uca4_methods <- c("PCA","cPCA++","Angry","Disgust","Split","Stack")

uca4_plot <- do.call(rbind,
lapply(seq_along(uca4_vectors), function(V) {
  res <- data.table(sus_dis_ans %*% uca4_vectors[[V]][,1:2], rep(c("SUS","DIS","ANS"), rep(35,3)),uca4_methods[V])
  setnames(res, c("PC1", "PC2", "Class","Method"))
  }))
uca4_plot[,Method := factor(Method, levels = c("PCA","cPCA++","Angry", "Disgust","Split","Stack"))]
ggplot(data = uca4_plot, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(title = "Female Angry, Disgust, Surprise ", x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("example/Faces/F_Angry_Disgust_Surprise_projected.png",  width = 11, height = 8, units = "in")


