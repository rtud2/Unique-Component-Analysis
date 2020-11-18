rm(list = ls())
Sys.setenv(MKL_DEBUG_CPU_TYPE=5)#, MKL_VERBOSE=1) #intel MKL settings. latter to hack amd processors
#setwd("/mnt/d/Residual-Dimension-Reduction")
setwd("D:/Residual-Dimension-Reduction/")
#install.packages("uca_0.13.zip", repos = NULL, type="source")
#install.packages("uca_0.13.tar.gz",type="source", repos = NULL) 
libraries <- c("data.table", "MASS", "ggplot2", "splines","gridExtra", "uca", "imager", "future", "RSpectra", "microbenchmark", "Rfast")
invisible(lapply(libraries, library, character.only = T))
plan(multiprocess)

## Helper functions 
source("functions/helper_functions.R")
kdef <- "faces/KDEF_and_AKDEF/KDEF/"
kdef_A <- "faces/KDEF_and_AKDEF/A_kdef_files/"

# Practice Emotions
sus_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"SUS.JPG")))
sas_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"SAS.JPG")))
nes_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"NES.JPG")))
dis_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"DIS.JPG")))
ans_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"ANS.JPG")))
has_pract_files <- paste0(kdef_A, sapply(list.files(kdef_A), function(zz) paste0(zz,"/",zz,"HAS.JPG")))

practice_list <- list(sus_pract_files, has_pract_files,sas_pract_files,nes_pract_files,dis_pract_files,ans_pract_files)

# practice_male_list <- lapply(practice_list, function(zz){grep(pattern = "AM", x = zz, value = T)})
# male_emotion_list %<-% lapply(practice_male_list, function(x) readImage(x, 50, 68))
# male_emotion_list <- lapply(male_emotion_list, transpose)
practice_female_list <- lapply(practice_list, function(zz){grep(pattern = "AF", x = zz, value = T)})
female_emotion_list %<-% lapply(practice_female_list, function(x) readImage(x, 50, 68))
female_emotion_list <- lapply(female_emotion_list, transpose)

# Final Emotions
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

# Stack Background Variability vs. BG1 and BG2 
# 1. sus & dis, 2. sus & sas, 3. nes & dis appear to be the most "different" from the stack
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
head(background_var[order(diff, decreasing = T)])

# bg1 bg2 Stack_var  bg1_var  bg2_var      diff
# 1: sus dis  786.1118 896.8199 734.6812 162.13865
# 2: sus sas  829.5219 911.0035 780.1336 130.86992
# 3: nes dis  770.3364 863.9260 754.1740 109.75200
# 4: sus has  841.1609 910.9744 826.0987  84.87568
# 5: sas nes  811.0090 793.3120 878.0899  84.77798
# 6: sus ans  817.5896 890.6497 807.1156  83.53412

# INNER PRODUCT Betwwen Stack, BG1 and BG2
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
head(max_inner[order(max)])
# bg1 bg2       max
# 1: sus ans 0.9691885
# 2: nes ans 0.9727309
# 3: nes dis 0.9791950
# 4: sus dis 0.9795358
# 5: ans dis 0.9812969
# 6: has ans 0.9814274
# 
background_inner[order(bg1_bg2)]
# bg1 bg2   Stack_1   Stack_2   bg1_bg2
# 1: sus ans 0.9691885 0.9665842 0.8800309
# 2: nes ans 0.9659761 0.9727309 0.8864153
# 3: sus dis 0.9795358 0.9586367 0.8959284
# 4: nes dis 0.9791950 0.9742804 0.9166890
# 5: sus sas 0.9864356 0.9676066 0.9167276
# 6: has ans 0.9766456 0.9814274 0.9220293

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


getFaceProjection=function(idx1, idx2){
  emo_names <- c("SUS","HAS", "SAS", "NES", "DIS","ANS")
  # function to cat class names onto output...
  bg_idx <- c(idx1,idx2)
  class_names = lapply(emo_names[-bg_idx], function(others){rep(c(others,emo_names[bg_idx]), each = 35)})  
  
  uca_bg1_bg2 <- lapply(final_female_emotion_list[-c(idx1,idx2)], function(targ){
    tmp_targ <- scale(do.call(rbind, c(list(targ), final_female_emotion_list[c(idx1,idx2)] )))
    tmp_bg1 <- scale(female_emotion_list[[idx1]])
    tmp_bg2 <- scale(female_emotion_list[[idx2]])
    tmp_stack <- scale(do.call(rbind, female_emotion_list[c(idx1,idx2)]))
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

#nes-dis
# order: sus, has, sas, ans
#background: female_emotion_list[c(4,5)] (NES,DIS)
uca_nes_dis <- getFaceProjection(4,5)
#sus-dis
# order: has, sas, dis, ans
#background: female_emotion_list[c(1,5)] (SUS,DIS)
uca_sus_dis <- getFaceProjection(1,5)

#sus-sas
# order: has, nes, dis, ans
#background: female_emotion_list[c(1,3)] (SUS,SAS)
uca_sus_sas <- getFaceProjection(1,3)

#ans-dis
# order: sus, has, sas, nes
#background: female_emotion_list[c(5,6)] (DIS,ANS)
uca_ans_dis <- getFaceProjection(5,6)

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
#ans_dis
#list order: sus, has, sad, nes
#within list order: PCA, cPCA++, ANS, DIS, Split, Stack
ans_dis_emo_list <- c("Surprise","Happy","Sad","Neutral")
uca_ans_dis_plot <- uca_ans_dis[[1]]

ans_dis_plot <- lapply(seq_along(uca_ans_dis_plot), function(idx){
  do.call(rbind, lapply(uca_ans_dis_plot[[idx]][,unique(method)], function(inner){
    toEigenfaces(data.matrix(uca_ans_dis_plot[[idx]][method == inner, .SD, .SDcols = -"method"]), 50, 68, inner)}
    ))
  } )
#plot eigenfaces
lapply(seq_along(ans_dis_emo_list), function(zz){
  plotEigenfaces2(ans_dis_plot[[zz]],"")   
  ggsave(paste0("example/Faces/Angry_Disgust/F_Angry_Disgust_",ans_dis_emo_list[zz],".png"),  width = 8, height = 11, units = "in")
  })
#plot projected faces
lapply(seq_along(uca_ans_dis[[2]]), function(idx){
  tmp_dat <- uca_ans_dis[[2]][[idx]]
  tmp_dat[,Method := factor(Method, levels = c("PCA","cPCA++","ANS","DIS","Split","Stack"))]
  ggplot(data = tmp_dat, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(#title = paste0("Female Surprise, Angry, ", ans_dis_emo_list[idx]),
    x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)
        )
ggsave(paste0("example/Faces/Angry_Disgust/F_Angry_Disgust_",ans_dis_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
    }        
  )

#sus_sas
#list order: has, nes, dis, ans
#within list order: PCA, cPCA++, SUS, SAS, Split, Stack
sus_sas_emo_list <- c("Happy","Neutral","Disgust","Angry")
uca_sus_sas_plot <- uca_sus_sas[[1]]

sus_sas_plot <- lapply(seq_along(uca_sus_sas_plot), function(idx){
  do.call(rbind, lapply(uca_sus_sas_plot[[idx]][,unique(method)], function(inner){
    toEigenfaces(data.matrix(uca_sus_sas_plot[[idx]][method == inner, .SD, .SDcols = -"method"]), 50, 68, inner)}
    ))
  } )
#plot eigenfaces
lapply(seq_along(sus_sas_emo_list), function(zz){
  plotEigenfaces2(sus_sas_plot[[zz]],"")   
  ggsave(paste0("example/Faces/Surprise_Sad/F_Surprise_Sad_",sus_sas_emo_list[zz],".png"),  width = 8, height = 11, units = "in")
  })
#plot projected faces
lapply(seq_along(uca_sus_sas[[2]]), function(idx){
  tmp_dat <- uca_sus_sas[[2]][[idx]]
  tmp_dat[,Method := factor(Method, levels = c("PCA","cPCA++","SUS","SAS","Split","Stack"))]
  ggplot(data = tmp_dat, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(#title = paste0("Female Surprise, Angry, ", sus_sas_emo_list[idx]),
    x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)
        )
ggsave(paste0("example/Faces/Surprise_Sad/F_Surprise_Sad_",sus_sas_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
    }        
  )

#sus_dis
#list order: has, sas, nes, ans
#within list order: PCA, cPCA++, SUS, DIS, Split, Stack
sus_dis_emo_list <- c("Happy","Sad","Neutral","Angry")
uca_sus_dis_plot <- uca_sus_dis[[1]]

sus_dis_plot <- lapply(seq_along(uca_sus_dis_plot), function(idx){
  do.call(rbind, lapply(uca_sus_dis_plot[[idx]][,unique(method)], function(inner){
    toEigenfaces(data.matrix(uca_sus_dis_plot[[idx]][method == inner, .SD, .SDcols = -"method"]), 50, 68, inner)}
    ))
  } )
#plot eigenfaces
lapply(seq_along(sus_dis_emo_list), function(zz){
  plotEigenfaces2(sus_dis_plot[[zz]],"")   
  ggsave(paste0("example/Faces/Surprise_Disgust/F_Surprise_Disgust_",sus_dis_emo_list[zz],".png"),  width = 8, height = 11, units = "in")
  })
#plot projected faces
lapply(seq_along(uca_sus_dis[[2]]), function(idx){
  tmp_dat <- uca_sus_dis[[2]][[idx]]
  tmp_dat[,Method := factor(Method, levels = c("PCA","cPCA++","SUS","DIS","Split","Stack"))]
  ggplot(data = tmp_dat, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(#title = paste0("Female Surprise, Angry, ", sus_dis_emo_list[idx]),
    x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)
        )
ggsave(paste0("example/Faces/Surprise_Disgust/F_Surprise_Disgust_",sus_dis_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
    }        
  )

#nes_dis
#list order: has, sas, nes, ans
#within list order: PCA, cPCA++, SUS, DIS, Split, Stack
sus_dis_emo_list <- c("Happy","Sad","Neutral","Angry")
uca_sus_dis_plot <- uca_sus_dis[[1]]

sus_dis_plot <- lapply(seq_along(uca_sus_dis_plot), function(idx){
  do.call(rbind, lapply(uca_sus_dis_plot[[idx]][,unique(method)], function(inner){
    toEigenfaces(data.matrix(uca_sus_dis_plot[[idx]][method == inner, .SD, .SDcols = -"method"]), 50, 68, inner)}
    ))
  } )
#plot eigenfaces
lapply(seq_along(sus_dis_emo_list), function(zz){
  plotEigenfaces2(sus_dis_plot[[zz]],"")   
  ggsave(paste0("example/Faces/Surprise_Disgust/F_Surprise_Disgust_",sus_dis_emo_list[zz],".png"),  width = 8, height = 11, units = "in")
  })
#plot projected faces
lapply(seq_along(uca_sus_dis[[2]]), function(idx){
  tmp_dat <- uca_sus_dis[[2]][[idx]]
  tmp_dat[,Method := factor(Method, levels = c("PCA","cPCA++","SUS","DIS","Split","Stack"))]
  ggplot(data = tmp_dat, aes(x = PC1, y = PC2, color = Class))+
  geom_point()+
  facet_wrap(~Method, scale = "free")+
  labs(#title = paste0("Female Surprise, Angry, ", sus_dis_emo_list[idx]),
    x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)
        )
ggsave(paste0("example/Faces/Surprise_Disgust/F_Surprise_Disgust_",sus_dis_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
    }        
  )

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
  labs(#title = paste0("Female Surprise, Angry, ", sus_ans_emo_list[idx]),
    x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)
        )
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
ggsave(paste0("example/Faces/Neutral_Angry/F_Neutral_Angry_",nes_ans_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
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
ggsave(paste0("example/Faces/Neutral_Disgust/F_Neutral_Disgust_",nes_dis_emo_list[idx],"_projected.png"),  width = 11, height = 8, units = "in")
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

