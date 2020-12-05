rm(list = ls())
library(ggplot2)
library(tidyverse)
library(geigen)
library(RSpectra)
library(splines)
library(uca)
library(latex2exp)
library(data.table)
source("../../functions/helper_functions.R")
## source("functions/rPCA.R")
## ##############################################################
## mouse data
## ##############################################################
mouse = read_csv("Data_Cortex_Nuclear.csv")
mouse[is.na(mouse)] = 0
## remove duplicate rows and columns
mouse = distinct(mouse) %>% select(-which(duplicated(t(mouse))))

mouse_dt <- data.table(mouse)

bg = filter(mouse, Behavior == "C/S", Treatment %in% c("Saline")) %>%
  select(-c(MouseID, class, Treatment, Behavior))
classes = as.factor(tg$Genotype)
B = scale(filter(bg, Genotype == "Control") %>% select(-Genotype))

# target: Saline Context-Shock
cs_gentype = "Ts65Dn"
#cs_gentype = "Control"
target_data <- data.matrix(mouse_dt[Treatment == "Saline" & Behavior == "C/S",lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
target_data_labs <-mouse_dt[Treatment == "Saline" & Behavior == "C/S", .SD, .SDcols = c("Genotype")]

bg1 <-data.matrix(mouse_dt[Treatment == "Memantine" & Behavior == "S/C" & Genotype == cs_gentype, lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
bg2 <- data.matrix(mouse_dt[Treatment == "Memantine" & Behavior == "C/S" & Genotype == cs_gentype, lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])
bg3 <- data.matrix(mouse_dt[Treatment == "Saline" & Behavior == "S/C" & Genotype ==cs_gentype, lapply(.SD, scale), .SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])

bg12 <- data.matrix(mouse_dt[Treatment == "Memantine"  & Genotype == cs_gentype, lapply(.SD, scale),.SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])

bg123 <- rbind(mouse_dt[Treatment == "Memantine"  & Genotype == cs_gentype],
               mouse_dt[Treatment == "Saline" & Behavior == "S/C" & Genotype == cs_gentype])
bg123 <- data.matrix(bg123[, lapply(.SD, scale),.SDcols = -c("MouseID","Genotype","Treatment","Behavior","class")])

target_cov <- cov(target_data)
bg1_cov <- cov(bg1)
bg2_cov <- cov(bg2)
bg3_cov <- cov(bg3)

bg12_stack_cov <- cov(bg12)
bg_stack_cov <- cov(bg123)

# show that cpca results cannot do better than autochosen results by UCA
cpca_bg1_0p5 <- eigs_sym(target_cov - 0.5*bg1_cov, which = "LA", k = 2)$vectors
cpca_bg1_1 <- eigs_sym(target_cov - 1*bg1_cov, which = "LA", k = 2)$vectors
cpca_bg1_5 <- eigs_sym(target_cov - 5*bg1_cov, which = "LA", k = 2)$vectors
cpca_bg1_10 <- eigs_sym(target_cov - 10*bg1_cov, which = "LA", k = 2)$vectors

cpca_bg2_0p5 <- eigs_sym(target_cov - 0.5*bg2_cov, which = "LA", k = 2)$vectors
cpca_bg2_1 <- eigs_sym(target_cov - 1*bg2_cov, which = "LA", k = 2)$vectors
cpca_bg2_5 <- eigs_sym(target_cov - 5*bg2_cov, which = "LA", k = 2)$vectors
cpca_bg2_10 <- eigs_sym(target_cov - 10*bg2_cov, which = "LA", k = 2)$vectors

cpca_bg3_0p5 <- eigs_sym(target_cov - 0.5*bg3_cov, which = "LA", k = 2)$vectors
cpca_bg3_1 <- eigs_sym(target_cov - 1*bg3_cov, which = "LA", k = 2)$vectors
cpca_bg3_5 <- eigs_sym(target_cov - 5*bg3_cov, which = "LA", k = 2)$vectors
cpca_bg3_10 <- eigs_sym(target_cov - 10*bg3_cov, which = "LA", k = 2)$vectors

cpca_stack_0p5 <- eigs_sym(target_cov - 0.5*bg_stack_cov, which = "LA", k = 2)$vectors
cpca_stack_1 <- eigs_sym(target_cov - 1*bg_stack_cov, which = "LA", k = 2)$vectors
cpca_stack_5 <- eigs_sym(target_cov - 5*bg_stack_cov, which = "LA", k = 2)$vectors
cpca_stack_10 <- eigs_sym(target_cov - 10*bg_stack_cov, which = "LA", k = 2)$vectors

# rpc_bg1 <- rPCA(target_data, bg1, n_components = 2, standardize = F)
# rpc_bg2 <- rPCA(target_data, bg2, n_components = 2, standardize = F)
# rpc_bg3 <- rPCA(target_data, bg3, n_components = 2, standardize = F)
# rpc_bg123 <- rPCA(target_data, bg123, n_components = 2, standardize = F)

cpcapp_bg1 <- eigen(solve(bg1_cov) %*% target_cov, 2)$vector[,1:2]
cpcapp_bg2 <- eigen(solve(bg2_cov) %*% target_cov, 2)$vector[,1:2]
cpcapp_bg3 <- eigen(solve(bg3_cov) %*% target_cov, 2)$vector[,1:2]
cpcapp_bg123 <- eigen(solve(bg_stack_cov) %*% target_cov, 2)$vector[,1:2]

plot_cpca_res <- rbind(
  data.table(target_data %*% cpca_bg1_0p5, "0.5", paste0("Mem-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg1_1, "1", paste0("Mem-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg1_5, "5", paste0("Mem-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg1_10, "10", paste0("Mem-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg2_0p5, "0.5", paste0("Mem-C/S-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg2_1, "1", paste0("Mem-C/S-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg2_5, "5", paste0("Mem-C/S-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg2_10, "10", paste0("Mem-C/S-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg3_0p5, "0.5", paste0("Saline-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg3_1, "1", paste0("Saline-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg3_5, "5", paste0("Saline-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_bg3_10, "10", paste0("Saline-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpca_stack_0p5, "0.5", paste0("Stack All (",cs_gentype,")"), target_data_labs),
  data.table(target_data %*% cpca_stack_1, "1", paste0("Stack All (",cs_gentype,")"), target_data_labs),
  data.table(target_data %*% cpca_stack_5, "5", paste0("Stack All (",cs_gentype,")"), target_data_labs),
  data.table(target_data %*% cpca_stack_10, "10", paste0("Stack All (",cs_gentype,")"), target_data_labs),
#  data.table(rpc_bg1, "rPCA", paste0("Mem-S/C-",cs_gentype), target_data_labs),
#  data.table(rpc_bg2, "rPCA",paste0("Mem-C/S-",cs_gentype), target_data_labs),
#  data.table(rpc_bg3, "rPCA", paste0("Saline-S/C-",cs_gentype), target_data_labs),
#  data.table(rpc_bg123, "rPCA", paste0("Stack All (",cs_gentype,")"), target_data_labs),
  data.table(target_data %*% cpcapp_bg1, "cPCA++", paste0("Mem-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpcapp_bg2, "cPCA++",paste0("Mem-C/S-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpcapp_bg3, "cPCA++", paste0("Saline-S/C-",cs_gentype), target_data_labs),
  data.table(target_data %*% cpcapp_bg123, "cPCA++", paste0("Stack All (",cs_gentype,")"), target_data_labs)
  )

setnames(plot_cpca_res, c("cPC1","cPC2","Parameter","Method","Label"))
plot_cpca_res[,Parameter := factor(Parameter,
                                   levels = c("cPCA++","rPCA","0.5", "1", "5", "10"))]
other_res <- ggplot(data = plot_cpca_res)+
  geom_point(aes(x = cPC1, y = cPC2, color = Label), alpha = 0.5)+
  labs(x = "Component 1", y = "Component 2")+
  facet_grid(Method~Parameter, scale = "free")+
  theme_bw(base_size = 13)+
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13), 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)
  )

ggsave(paste0("Mouse_stack_cpc_rpc",cs_gentype,".png"), other_res, width = 11, height = 8, units = "in")

pc <- eigs_sym(target_cov, which = "LA", k = 2)$vectors
uca_bg1 <- uca(A = target_cov, B = bg1_cov)$vectors
uca_bg2 <- uca(A = target_cov, B = bg2_cov)$vectors
uca_bg3 <- uca(A = target_cov, B = bg3_cov)$vectors
uca_stack <- uca(A = target_cov, B = bg_stack_cov, method = "cov")$vectors
uca_split <- uca(A = target_cov, B = list(bg1_cov, bg2_cov,bg3_cov), method = "cov")$vectors
cpca_pp <- eigen(solve(bg_stack_cov) %*% target_cov)$vectors[,1:2]

plot_memantine_saline <- rbind(
  data.table(target_data %*% pc, paste0("PCA"), target_data_labs),
  data.table(target_data %*% uca_bg1, paste0("UCA (Memantine-S/C-",cs_gentype,")"), target_data_labs),
  data.table(target_data %*% uca_bg2, paste0("UCA (Memantine-C/S-",cs_gentype,")"), target_data_labs),
  data.table(target_data %*% uca_bg3, paste0("UCA (Saline-S/C-",cs_gentype,")"), target_data_labs),
  data.table(target_data %*% cpca_pp, "cPCA++ (Stack All)", target_data_labs),
#  data.table(rpc_bg123, "rPCA (Stack All)", target_data_labs),
  data.table(target_data %*% uca_stack, "UCA (Stack All)", target_data_labs),
  data.table(target_data %*% uca_split, "UCA (Split All)", target_data_labs))
setnames(plot_memantine_saline, c("UC1","UC2","Method","Label"))
plot_memantine_saline[, Method := factor(Method, levels = unique(Method))]

mouse_split_stack <- ggplot(data = plot_memantine_saline)+
  geom_point(aes(x = UC1, y = UC2, color = Label), alpha = 0.5)+
  labs(x = "Component 1", y = "Component 2")+
  facet_wrap(~Method, scale = "free", nrow = 2)+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom")

ggsave(paste0("Mouse_split_stack_",cs_gentype,".png"), mouse_split_stack, width = 11, height = 8, units = "in")
