
readImage <- function(listFiles, size_x, size_y){
  loaded_images <- lapply(listFiles, load.image)
  grey_scaled <- lapply(loaded_images, function(zz) scale(as.data.frame(resize(grayscale(zz), size_x, size_y))[,3]))
  return(data.matrix(do.call(cbind,grey_scaled)))
}

# Turn Eigenvectors to eigenfaces
toEigenfaces <- function(dat, size_x = 58, size_y = 68, name){
  eigenfaces <- data.table(do.call(rbind,
                                   lapply(1:ncol(dat), function(zz){
                                     cbind(as.data.frame(as.cimg(dat[,zz], dim = c(size_x, size_y,1,1))), zz)
                                   })))
  setnames(eigenfaces,c("x","y","value","Eigenface"))
  eigenfaces[, alpha := name]
  eigenfaces[, value2 := 2*(value - min(value))/(max(value) - min(value)), by = Eigenface]
  return(eigenfaces)
}

plotEigenfaces = function(eigenfaces, title, nrow){
  ggplot(eigenfaces, aes(x,y))+
    geom_raster(aes(fill = value))+
    facet_wrap(~Eigenface, nrow = nrow)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())+
    scale_fill_gradient(low="black", high="white")+
    labs(title = title)
}

plotEigenfaces2 <- function(eigenfaces, title){
  ggplot(eigenfaces, aes(x,y))+
    geom_raster(aes(fill = value))+
    facet_grid(alpha~Eigenface)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())+
    scale_fill_gradient(low="black", high="white")+
    labs(title = title)+
    theme(strip.text = element_text(size = 20),
          axis.title = element_blank(),
          legend.position = "bottom",
          axis.ticks = element_blank(),
          axis.text =  element_blank())
}

readImageRotate <- function(listFiles, rotation, size_x, size_y){
  loaded_images <- lapply(listFiles, load.image)
  grey_scaled <- lapply(seq_along(loaded_images), function(zz) scale(as.data.frame(resize(grayscale(imrotate(loaded_images[[zz]], rotation[[zz]])), size_x, size_y))[,3]))
  return(data.matrix(do.call(cbind,grey_scaled)))
}

# make faces positively correlated with reference data
pos_corr <- function(dat, ref_dat){
  sapply(1:ncol(dat), function(zz){
    current_sign <- sign(cor(dat[,zz], ref_dat[,zz]))
    if(current_sign == -1){
      return(-dat[,zz])
    }else{
      return(dat[,zz])
    }
  })
}
