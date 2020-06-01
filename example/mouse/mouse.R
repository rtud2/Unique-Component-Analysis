rm(list = ls())

library(ggplot2)
library(tidyverse)
library(geigen)
library(RSpectra)
library(uca)
source("functions/helper_functions.R")

## ##############################################################
## mouse data
## ##############################################################
mouse = read_csv("../../../contrastive/experiments/datasets/Data_Cortex_Nuclear.csv")
mouse[is.na(mouse)] = 0

## remove duplicate rows and columns
mouse = distinct(mouse) %>% select(-which(duplicated(t(mouse))))

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
v = eigen(cov(X) - 1 * cov(B))$vectors[, 1:2]
cpca1 = X %*% v
v = eigen(cov(X) - 3 * cov(B))$vectors[, 1:2]
cpca3 = X %*% v
v = eigen(cov(X) - 5 * cov(B))$vectors[, 1:2]
cpca5 = X %*% v
v = eigen(cov(X) - 7 * cov(B))$vectors[, 1:2]
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
                Method = c(rep("UCA", n),
                           rep("PCA", n),
                           rep("cPCA alpha = 1", n),
                           rep("cPCA alpha = 3", n),
                           rep("cPCA alpha = 5", n),
                           rep("cPCA alpha = 7", n),
                           rep("rPCA", n),
                           rep("cPCA++", n)),
                Classes = rep(classes, 8))

mouse_comparison <- ggplot(df, aes(x = pc1, y = pc2, color = Classes)) +
  geom_point() +
  facet_wrap(~ Method, scale = "free", nrow = 2)+
  labs(x = "Component 1", y = "Component 2")
ggsave("example/Mouse_Data.png", mouse_comparison, width = 11, height = 8, units = "in")




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
  labs(x = "Component 1", y = "Component 2")

ggsave("example/Mouse_split_stack.png", mouse_split_stack, width = 11, height = 8, units = "in")
