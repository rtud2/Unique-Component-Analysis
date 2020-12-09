library("uca")
library("microbenchmark")
library("data.table")
library("ggplot2")

p = c(seq(from = 100, to = 1000, by = 100),seq(from = 2000, to = 10000, by = 1000))#, seq(from = 20E3, to = 50E3, by = 10E3))

res <- list()

for( i in seq_along(p)){
  

timing_data <- microbenchmark(data = uca(target, bg, method = "data", nv = 2),
                              cov= uca(cov_target, cov_bg, method = "cov", nv = 2),
                              unit = "s",
                              times = 25,
                              setup = {
  bg <- matrix(rnorm(100*p[i], mean = 0, sd = 1),ncol = p[i])
  target <- matrix(rnorm(100*p[i]),ncol = p[i])+bg
  cov_target <- cov(target)
  cov_bg <- cov(bg)
  
                              }
)
size <- c("data" = object.size(target), "cov" = object.size(cov_target))

tmp_res <- data.table(do.call(cbind, timing_data))
setnames(tmp_res, c("method", "time"))

tmp_res[, `:=`(method = timing_data$expr,
               size = size[method],
               dim  = p[i])]

res[[i]] <- tmp_res

res_dt <- do.call(rbind, res)
#saveRDS(res, "example/res.rds")
fwrite(res_dt, "example/timing_res.csv")
res_dt <- fread("timing_res.csv")

current <- ggplot(data = res_dt, aes(x = as.factor(dim), y = (time/10e9), fill = as.factor(method)))+
                  geom_boxplot(alpha = 0.5)+
                  labs(x = "Dimension", y = "Time (Sec.)")+
                  theme_bw(base_size = 20)+
                  theme(legend.position = "bottom", axis.text.x = element_text(angle = -45))

ggsave(filename = paste0("example/current_",p[i],".png"),plot = current, width = 11, height = 8, units = "in")

message(paste("done with: ", p[i], "\n"))
}

res_dt[, method := factor(method, levels = c("cov","data"), labels = c("Eigendecomposition","Product SVD"))]

final <- ggplot(data = res_dt[dim >=1000], aes(x = as.factor(dim), y = (time/10e9), fill = method))+
                  geom_boxplot(alpha = 0.5)+
                  labs(x = "Dimension", y = "Time (Sec.)")+
                  theme_bw(base_size = 20)+
                  scale_fill_discrete(name = "Method")+
                  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0))

ggsave(filename = "final_perf.png",plot = final, width = 11, height = 8, units = "in")



final_log <- ggplot(data = res_dt[dim >=1000], aes(x = as.factor(dim), y = log(time/10e6), fill = method))+
                  geom_boxplot(alpha = 0.5)+
                  labs(x = "Dimension", y = "log(Time) (milliseconds)")+
                  theme_bw(base_size = 20)+
                  scale_fill_discrete(name = "Method")+
                  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0))

ggsave(filename = "final_perf_log.png",plot = final_log, width = 11, height = 8, units = "in")
