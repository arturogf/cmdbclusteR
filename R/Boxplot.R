
stayboxplotfile <- file.path(directory, "data/output", paste(nombre_generado, "LOS-boxplot", ".pdf", sep=""))

library("ggplot2")

pdf(stayboxplotfile, paper = "a4r")

#plot the boxplot
a<-ggplot(salida,
       aes(
         group = salida[,pos_cluster],
         x = salida[,pos_cluster],
         y = salida[, pos_stay]
       ))  +
  coord_flip() +
  geom_boxplot(
    #aes(fill = diagnostico1),
    #fill="white",
    colour = "black",
    outlier.colour = "red",
    outlier.shape = 2
  ) +
  stat_boxplot(geom = 'errorbar') +
  geom_jitter(width = 0.2, height=0.2) +
  #geom_point(position=position_dodge(width=0.75),aes(group=salida$cluster))+
  labs(title = "Boxplot: LOS per cluster") +
  ylab("Length of Stay(LOS)") +
  scale_x_continuous(name = "Num. cluster",
                     breaks = seq(1:length(unique(salida[,pos_cluster]))),
                     labels = labls)
print(a)
dev.off()