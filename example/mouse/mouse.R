rm(list = ls())

library(ggplot2)
library(tidyverse)
library(geigen)
library(RSpectra)
library(uca)
library(latex2exp)
library(data.table)
source("functions/helper_functions.R")

## ##############################################################
## mouse data
## ##############################################################
mouse = read_csv("../contrastive/experiments/datasets/Data_Cortex_Nuclear.csv")
mouse[is.na(mouse)] = 0

## remove duplicate rows and columns
mouse = distinct(mouse) %>% select(-which(duplicated(t(mouse))))

mouse_dt <- data.table(mouse)

tg = filter(mouse, Behavior == "S/C", Treatment %in% c("Saline")) %>%
  select(-c(MouseID, class, Treatment, Behavior))
bg = filter(mouse, Behavior == "C/S", Treatment %in% c("Saline")) %>%
  select(-c(MouseID, class, Treatment, Behavior))

classes = as.factor(tg$Genotype)
X = scale(select(tg, -Genotype))
B = scale(filter(bg, Genotype == "Control") %>% select(-Genotype))
B2 = scale(filter(bg, Genotype == "Ts65Dn") %>% select(-Genotype))

## ##############################################################
## uca
## ##############################################################
covX = cov(X)
covB = cov(B)
covB2 = cov(B2)

uca_res <- X %*% uca(X, B, nv = 2)$vectors
## ##############################################################
## standard pca
## ##############################################################
v = svd(X)$v[, 1:2]
pca = X %*% v

## ##############################################################
## rpca
## ##############################################################
v_B = svd(B)$v[, 1:2]
v = svd(X %*% (diag(ncol(B)) - tcrossprod(v_B, v_B)))$v[, 1:2]
rpca = X %*% v

## ##############################################################
## cpca
## ##############################################################
v = eigen(cov(X) - 0.5 * cov(B))$vectors[, 1:2]
cpca1 = X %*% v
v = eigen(cov(X) - 1 * cov(B))$vectors[, 1:2]
cpca3 = X %*% v
v = eigen(cov(X) - 5 * cov(B))$vectors[, 1:2]
cpca5 = X %*% v
v = eigen(cov(X) - 10 * cov(B))$vectors[, 1:2]
cpca7 = X %*% v


## ##############################################################
## cpca++
## ##############################################################
v = eigen(solve(covB) %*% covX)$vectors[, 1:2]
cpcapp = X %*% v

crossprod(v) # non-orthogonal

## ##############################################################
## visualize
## ##############################################################
pca = pos_corr(pca, uca_res)
cpca1 = pos_corr(cpca1, uca_res)
cpca3 = pos_corr(cpca3, uca_res)
cpca5 = pos_corr(cpca5, uca_res)
cpca7 = pos_corr(cpca7, uca_res)
rpca = pos_corr(rpca, uca_res)
cpcapp = pos_corr(cpcapp, uca_res)

n = nrow(X)
df = data.frame(pc1 = c(uca_res[,1], pca[,1], cpca1[,1], cpca3[,1], cpca5[,1], cpca7[,1], rpca[,1], cpcapp[,1]),
                pc2 = c(uca_res[,2], pca[,2], cpca1[,2], cpca3[,2], cpca5[,2], cpca7[,2], rpca[,2], cpcapp[,2]),
                Method = factor(c(rep("UCA", n),
                           rep("PCA", n),
                           rep("cPCA alpha = 0.5", n),
                           rep("cPCA alpha = 1", n),
                           rep("cPCA alpha = 5", n),
                           rep("cPCA alpha = 10", n),
                           rep("rPCA", n),
                           rep("cPCApp", n)),
                           levels = c("cPCA alpha = 0.5", "cPCA alpha = 1", "cPCA alpha = 5", "cPCA alpha = 10", "PCA", "rPCA", "cPCApp","UCA"),
                           labels = c(TeX("cPCA $\\alpha =0.5$"), TeX("cPCA $\\alpha =1$"),TeX("cPCA $\\alpha =5$"), TeX("cPCA $\\alpha =10$") , "PCA", "rPCA", TeX("cPCA$++$"),"UCA") ),
                           Classes = rep(classes, 8))

# new_lab <- as_labeller(c(`cPCA alpha = 0.5` = bquote("cPCA"~alpha==~"0.5"),
#                          `cPCA alpha = 1` = bquote("cPCA"~alpha==~"1"),
#                          `cPCA alpha = 5` = bquote("cPCA"~alpha==~"5"),
#                          `cPCA alpha = 10` = bquote("cPCA"~alpha==~"10"),
#                          `PCA` = "PCA", `rPCA`= "rPCA", `cPCA++` = "cPCA++",`UCA` = "UCA"), label_bquote)

mouse_comparison <- ggplot(df, aes(x = pc1, y = pc2, color = Classes)) +
  geom_point() +
  facet_wrap(~ Method, scale = "free", nrow = 2, labeller = label_parsed)+
  labs(x = "Component 1", y = "Component 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)
        )

ggsave("Mouse_Data.png", mouse_comparison, width = 11, height = 8, units = "in")




mouse_stack <- uca(X, rbind(B, B2), nv = 2)
mouse_split <- uca(X, list(B, B2), nv = 2)

uca_stack <- X %*% mouse_stack$vectors
uca_split <- X %*% mouse_split$vectors

uca_stack = pos_corr(uca_stack, uca_split)

df2 <- data.frame(pc1 = c(uca_stack[,1], uca_split[,1]),
                  pc2 = c(uca_stack[,2], uca_split[,2]),
                  Method = c(rep("Stack", n),
                             rep("Split", n)),
                  Classes = rep(classes, 6))

mouse_split_stack <- ggplot(df2, aes(x = pc1, y = pc2, color = Classes)) +
  geom_point() +
  facet_wrap(~ Method, scale = "free")+
  labs(x = "Component 1", y = "Component 2")+
  theme_bw()

ggsave("Mouse_split_stack.png", mouse_split_stack, width = 11, height = 8, units = "in")


b_control_all <- mouse_dt[Behavior == "C/S" , lapply(.SD, scale), .SDcols = -c("MouseID", "class", "Treatment", "Behavior","Genotype")]
b_control_saline <- mouse_dt[Behavior == "C/S" & Treatment == "Saline", lapply(.SD, scale), .SDcols = -c("MouseID", "class", "Treatment", "Behavior","Genotype")]
b_control_memantine <- mouse_dt[Behavior == "C/S" & Treatment == "Memantine", lapply(.SD, scale), .SDcols = -c("MouseID", "class", "Treatment", "Behavior", "Genotype")]

t_shock_all <- mouse_dt[Behavior == "S/C", lapply(.SD, scale), .SDcols = -c("MouseID", "class", "Treatment", "Behavior", "Genotype")]
t_shock_saline <- mouse_dt[Behavior == "S/C" & Treatment == "Saline", lapply(.SD, scale), .SDcols = -c("MouseID", "class", "Treatment", "Behavior", "Genotype")]
t_shock_memantine <- mouse_dt[Behavior == "S/C" & Treatment == "Memantine", lapply(.SD, scale), .SDcols = -c("MouseID", "class", "Treatment", "Behavior", "Genotype")]

t_shock_all_labs <- mouse_dt[Behavior == "S/C", Genotype]
t_shock_memantine_labs <- mouse_dt[Behavior == "S/C" & Treatment == "Memantine", Genotype]
t_shock_saline_labs <- mouse_dt[Behavior == "S/C" & Treatment == "Saline", Genotype]
uca_memantine <-data.matrix(t_shock_memantine) %*% uca(data.matrix(t_shock_memantine), data.matrix(b_control_memantine))$vectors
uca_memantine_2 <-data.matrix(t_shock_memantine) %*% uca(data.matrix(t_shock_memantine), list(data.matrix(b_control_memantine), data.matrix(b_control_saline)))$vectors

cov_shock_saline <- cov(t_shock_saline)
cov_shock_memantine <- cov(t_shock_memantine)
cov_shock_all <- cov(t_shock_all)
cov_control_all <- cov(b_control_all)
cov_control_saline <- cov(b_control_saline)
cov_control_memantine <- cov(b_control_memantine)

cov_shock_m_control <- cov(rbind(t_shock_memantine, b_control_all))

#shock saline ~ vs. control saline and mematine
uca_saline_split <- data.matrix(t_shock_saline) %*% uca(cov_shock_saline, list(cov_control_saline, cov_control_memantine, cov_shock_memantine), method = "cov")$vectors 
uca_saline_stack <- pos_corr(data.matrix(t_shock_saline) %*% uca(cov_shock_saline,cov_shock_m_control, method = "cov")$vectors, uca_saline_split)
split_cpca0p1 <- pos_corr(data.matrix(t_shock_saline) %*% eigen(cov_shock_saline - 0.1*cov_shock_m_control)$vectors[,1:2], uca_saline_split)
split_cpca1 <- pos_corr(data.matrix(t_shock_saline) %*% eigen(cov_shock_saline - cov_shock_m_control)$vectors[,1:2], uca_saline_split)
split_cpcapp <- pos_corr(data.matrix(t_shock_saline) %*% eigen(solve(cov_shock_m_control)%*% cov_shock_saline )$vectors[,1:2], uca_saline_split)
split_cpca10 <- pos_corr(data.matrix(t_shock_saline) %*% eigen(cov_shock_saline - 10*cov_shock_m_control)$vectors[,1:2], uca_saline_split)

plot_saline <- rbind(
      data.table(uca_saline_split, "UCA Split", labs = t_shock_saline_labs),
      data.table(split_cpca0p1, "cPCA Stack 0.1", labs = t_shock_saline_labs),
      data.table(split_cpca1, "cPCA Stack 1", labs = t_shock_saline_labs),
      data.table(split_cpcapp, "cPCA++ Stack", labs = t_shock_saline_labs),
      data.table(split_cpca10, "cPCA Stack 10", labs = t_shock_saline_labs),
      data.table(uca_saline_stack, "UCA Stack", labs = t_shock_saline_labs))

setnames(plot_saline, c("cPC1", "cPC2", "Method", "Label"))

ggplot(data = plot_saline)+
geom_point(aes(cPC1, cPC2, color = Label), alpha = 0.5)+
  facet_wrap(~Method, scale = "free")+
  theme_bw()




# shock mematine ~ vs. control saline and mematine 
uca_all <- data.matrix(t_shock_memantine) %*% uca(cov_shock_memantine, list(cov_control_saline, cov_control_memantine, cov_shock_saline), method = "cov")$vectors 
uca_all_stack <- pos_corr(data.matrix(t_shock_memantine) %*% uca(cov_shock_memantine,cov_control_all, method = "cov")$vectors, uca_all)
all_cpca0p1 <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(cov_shock_memantine - 0.1*cov_control_all)$vectors[,1:2], uca_all)
all_cpca1 <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(cov_shock_memantine - cov_control_all)$vectors[,1:2], uca_all)
all_cpcapp <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(solve(cov_control_all)%*% cov_shock_memantine )$vectors[,1:2], uca_all)
all_cpca10 <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(cov_shock_memantine - 10*cov_control_all)$vectors[,1:2], uca_all)

plot_memantine <- rbind(
      #data.table(uca_memantine, "Control Memantine", labs = t_shock_memantine_labs),
      #data.table(uca_memantine_2, "Control Memantine & Saline", labs = t_shock_memantine_labs),
      data.table(uca_all, "UCA Split", labs = t_shock_memantine_labs),
      data.table(all_cpca0p1, "cPCA Stack 0.1", labs = t_shock_memantine_labs),
      data.table(all_cpca1, "cPCA Stack 1", labs = t_shock_memantine_labs),
      data.table(all_cpcapp, "cPCA++ Stack", labs = t_shock_memantine_labs),
      data.table(all_cpca10, "cPCA Stack 10", labs = t_shock_memantine_labs),
      data.table(uca_all_stack, "UCA Stack", labs = t_shock_memantine_labs))

setnames(plot_memantine, c("cPC1", "cPC2", "Method", "Label"))

ggplot(data = plot_memantine)+
geom_point(aes(cPC1, cPC2, color = Label), alpha = 0.5)+
  facet_wrap(~Method, scale = "free")+
  theme_bw()


#shock menantine vs shock saline and control menantine 
uca_all <- data.matrix(t_shock_memantine) %*% uca(cov_shock_memantine, list(covX, cov_control_memantine), method = "cov")$vectors 
uca_all_stack <- pos_corr(data.matrix(t_shock_memantine) %*% uca(cov_shock_memantine,covX)$vectors, uca_all)
all_cpca0p1 <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(cov_shock_memantine - 0.1*covX)$vectors[,1:2], uca_all)
all_cpca1 <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(cov_shock_memantine - covX)$vectors[,1:2], uca_all)
all_cpcapp <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(solve(covX)%*% cov_shock_memantine )$vectors[,1:2], uca_all)
all_cpca10 <- pos_corr(data.matrix(t_shock_memantine) %*% eigen(cov_shock_memantine - 10*covX)$vectors[,1:2], uca_all)

plot_memantine <- rbind(
      #data.table(uca_memantine, "Control Memantine", labs = t_shock_memantine_labs),
      #data.table(uca_memantine_2, "Control Memantine & Saline", labs = t_shock_memantine_labs),
      data.table(uca_all, "UCA Split", labs = t_shock_memantine_labs),
      data.table(all_cpca0p1, "cPCA Stack 0.1", labs = t_shock_memantine_labs),
      data.table(all_cpca1, "cPCA Stack 1", labs = t_shock_memantine_labs),
      data.table(all_cpcapp, "cPCA++ Stack", labs = t_shock_memantine_labs),
      data.table(all_cpca10, "cPCA Stack 10", labs = t_shock_memantine_labs),
      data.table(uca_all_stack, "UCA Stack", labs = t_shock_memantine_labs))

setnames(plot_memantine, c("cPC1", "cPC2", "Method", "Label"))

ggplot(data = plot_memantine)+
geom_point(aes(cPC1, cPC2, color = Label), alpha = 0.5)+
  facet_wrap(~Method, scale = "free")+
  theme_bw()

#target context-shock mice (both treatments), background shock-context mice (split vs stack)
target_sc_dat <- data.matrix(mouse_dt[Behavior == "S/C", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
target_sc <- cov(target_sc_dat)
bg_cs_stack <- cov(data.matrix(mouse_dt[Behavior == "C/S", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")]))
bg_cs_saline <- cov(data.matrix(mouse_dt[Behavior == "C/S" & Treatment == "Saline", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")]))
bg_cs_memantine <- cov(data.matrix(mouse_dt[Behavior == "C/S" & Treatment != "Saline", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")]))

target_sc_labs <- mouse_dt[Behavior == "S/C", .SD, .SDcols = "Genotype"]

sc_cs_stack <- uca(A = target_sc, B = bg_cs_stack, method = "cov")$vectors
sc_cs_split <- uca(A = target_sc, B = list(bg_cs_saline, bg_cs_memantine), method = "cov")$vectors

plot_cs_sc <- rbind(data.table(target_sc_dat %*% sc_cs_stack, "Stack", target_sc_labs),
                    data.table(target_sc_dat %*% sc_cs_split, "Split", target_sc_labs))
setnames(plot_cs_sc, c("UC1","UC2","Method","Label"))
ggplot(data = plot_cs_sc)+
  geom_point(aes(x = UC1, y = UC2, color = Label), alpha = 0.5)+
  facet_wrap(~Method, scale = "free")+
  theme_bw()

#target memantine (C/S and S/C), background saline (C/S and S/C).
target_memantine <- data.matrix(mouse_dt[Treatment == "Memantine", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
target_memantine_cov <- cov(target_memantine)
bg_saline_stack <- cov(data.matrix(mouse_dt[Treatment == "Saline", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")]))
bg_saline_sc <- cov(data.matrix(mouse_dt[Treatment == "Saline" & Behavior == "S/C", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")]))
bg_saline_cs <- cov(data.matrix(mouse_dt[Treatment == "Saline" & Behavior == "C/S", .SD, .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")]))

target_memantine_labs <- mouse_dt[Treatment == "Memantine", .SD, .SDcols = c("Genotype")]

memantine_stack <- uca(A = target_memantine_cov, B = bg_saline_stack, method = "cov")$vectors
memantine_split <- uca(A = target_memantine_cov, B = list(bg_saline_sc,bg_saline_cs), method = "cov")$vectors

plot_memantine_saline <- rbind(data.table(target_memantine %*% memantine_stack, "Stack", target_memantine_labs),
                    data.table(target_memantine %*% memantine_split, "Split", target_memantine_labs))
setnames(plot_memantine_saline, c("UC1","UC2","Method","Label"))
ggplot(data = plot_memantine_saline)+
  geom_point(aes(x = UC1, y = UC2, color = Label), alpha = 0.5)+
  facet_wrap(~Method, scale = "free")+
  theme_bw()

# target: Saline Shock-Context
target_data <- data.matrix(mouse_dt[Treatment == "Saline" & Behavior == "S/C",lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
target_data_labs <-mouse_dt[Treatment == "Saline" & Behavior == "S/C", .SD, .SDcols = c("Genotype")]

bg1 <-data.matrix(mouse_dt[Treatment == "Memantine" & Behavior == "S/C" & Genotype == "Control", lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
bg2 <- data.matrix(mouse_dt[Treatment == "Memantine" & Behavior == "C/S" & Genotype == "Control", lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
bg3 <- data.matrix(mouse_dt[Treatment == "Saline" & Behavior == "C/S" & Genotype == "Control", lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])

bg123 <- rbind(mouse_dt[Treatment == "Memantine"  & Genotype == "Control"],
               mouse_dt[Treatment == "Saline" & Behavior == "C/S" & Genotype == "Control"])
bg123 <- data.matrix(bg123[, lapply(.SD, scale),.SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])

target_cov <- cov(target_data)
bg1_cov <- cov(bg1)
bg2_cov <- cov(bg2)
bg3_cov <- cov(bg3)
#bg_stack_cov <- cov(bg23)
bg_stack_cov <- cov(bg123)

uca_bg1 <- uca(A = target_cov, B = bg1_cov)$vectors
uca_bg2 <- uca(A = target_cov, B = bg2_cov)$vectors
uca_bg3 <- uca(A = target_cov, B = bg3_cov)$vectors
uca_stack <- uca(A = target_cov, B = bg_stack_cov, method = "cov")$vectors
uca_split <- uca(A = target_cov, B = list(bg1_cov, bg2_cov,bg3_cov), method = "cov")$vectors

plot_memantine_saline <- rbind(
  data.table(target_data %*% uca_bg1, "BG1", target_data_labs),
  data.table(target_data %*% uca_bg2, "BG2", target_data_labs),
  data.table(target_data %*% uca_bg3, "BG3", target_data_labs),
  data.table(target_data %*% uca_stack, "Stack All", target_data_labs),
  data.table(target_data %*% uca_split, "Split All", target_data_labs))
setnames(plot_memantine_saline, c("UC1","UC2","Method","Label"))
ggplot(data = plot_memantine_saline)+
  geom_point(aes(x = UC1, y = UC2, color = Label), alpha = 0.5)+
  facet_wrap(~Method, scale = "free")+
  theme_bw()
  
