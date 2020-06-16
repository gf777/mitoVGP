library(dplyr)
library(ggplot2)
library(broom)
library(lawstat)
library("readxl")

setwd("~/Documents/VGP/com/papers/mitoVGP/submission_version/LM/")

mitoVGP <- read_excel("../Supplementary Tables/Supplementary Tables.xlsx", sheet = "ST1")

is.factor(mitoVGP$`Size selection kbp (Pacbio)`)
is.numeric(mitoVGP$`Total raw data (Gbp)`)

mitoVGP <- mutate(mitoVGP, `Size selection cutoff` = ifelse(`Size selection kbp (Pacbio)` > 20, "Above 20", "Below 20"))

mitoVGP$`Size selection cutoff` <- as.factor(mitoVGP$`Size selection cutoff`)
is.factor(mitoVGP$`Size selection cutoff`)

out1<-with(mitoVGP, table(`Tissue type (Pacbio)`, `Taxonomic group`))
out2<-with(mitoVGP, table(`Tissue type (Pacbio)`, Success))
out3<-with(mitoVGP, table(`Size selection kbp (Pacbio)`, Success))
out4<-with(mitoVGP, table(`Size selection cutoff`, Success))

write.table(out1, file='../Supplementary Tables/ST3.tsv', quote=FALSE, sep='\t')
write.table(out2, file='../Supplementary Tables/ST5.tsv', quote=FALSE, sep='\t')
write.table(out3, file='../Supplementary Tables/ST6a.tsv', quote=FALSE, sep='\t')
write.table(out4, file='../Supplementary Tables/ST6b.tsv', quote=FALSE, sep='\t')

mitoVGP <- mitoVGP %>%
  mutate(Success = ifelse(Success == "no",0,1))

png("../Extended Data Fig. 3/Extended Data Fig. 3.png", width = 2000, height = 750)

ggplot(mitoVGP, aes(x=reorder(`Tissue type (Pacbio)`, -`Available Pacbio mtDNA reads`),  y=`Available Pacbio mtDNA reads`)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=8,
               outlier.size=4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+
  ylab("Number of Pacbio reads")+
  theme(
    plot.margin = margin(t = 20, b = 0),
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
  )

dev.off()

plot(mitoVGP$Success, mitoVGP$`Available Pacbio mtDNA reads`)
plot(mitoVGP$`Taxonomic group`, mitoVGP$`Available Pacbio mtDNA reads`)

#Q <- quantile(mitoVGP$N.Pacbio.reads, probs=c(.25, .75), na.rm = FALSE)
#iqr <- IQR(mitoVGP$N.Pacbio.reads)
#mitoVGP_nooutliers<- subset(mitoVGP, mitoVGP$N.Pacbio.reads > (Q[1] - 1.5*iqr) & mitoVGP$N.Pacbio.reads < (Q[2]+1.5*iqr))

mod.mitoVGP <- glm(Success ~ `Available Pacbio mtDNA reads`, family="binomial", data=mitoVGP)
summary(mod.mitoVGP)
anova(mod.mitoVGP, test="Chi")

png("../Extended Data Fig. 2/Extended Data Fig. 2.png", width = 2000, height = 750)

mitoVGP %>%
  mutate(size = ifelse(`Available Pacbio mtDNA reads` < 1, "0", ">0")) %>%
  ggplot(aes(x = `Available Pacbio mtDNA reads`, y = Success, colour = size))  +
  stat_sum() +
  geom_smooth(method = "glm", formula = y ~ x, se = FALSE,
              method.args = list(family = binomial()), colour = "gray")+
  scale_x_continuous(trans='sqrt', breaks = c(5,100,200,300,400,500))+
  scale_color_manual(values = c('green', 'red'))+
  scale_y_continuous(name="Success", breaks = c(0, 1), labels = c("Failure", "Success"))+
  scale_size_continuous(range = c(5, 15))+
  xlab("Number of Pacbio reads")+
  theme(
    plot.margin = margin(t = 20, b = 0),
    axis.title.x = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.title.y = element_blank(),
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    legend.text=element_text(size=30),
    legend.title=element_blank()
  )+
  guides(colour = guide_legend(override.aes = list(size=10)))

dev.off()

#r2
nullmod <- glm(Success ~ 0, family="binomial", data=mitoVGP)
1-logLik(mod.mitoVGP)/logLik(nullmod)

#dist
hist(mitoVGP$`Available Pacbio mtDNA reads`)
hist(sqrt(mitoVGP$`Available Pacbio mtDNA reads`))

filtered_mitVGP<- mitoVGP %>% group_by(`Tissue type (Pacbio)`) %>% filter(n() > 3) %>%
                              group_by(`Taxonomic group`) %>% filter(n() > 3) %>%
                              group_by(`DNA extraction (Pacbio)`) %>% filter(n() > 3) %>%
                              group_by(`Library prep - fragmentation (Pacbio)`) %>% filter(n() > 3) %>%
                              group_by(`Library prep (Pacbio)`) %>% filter(n() > 3)
  

#N reads by tissue type, taxonomic group, size selection, DNA extraction, Library prep fragmentation, library prep
mod.mitoVGP <- lm(`Available Pacbio mtDNA reads` ~ `Tissue type (Pacbio)` + `Taxonomic group` + `Size selection kbp (Pacbio)` + `DNA extraction (Pacbio)` + `Library prep - fragmentation (Pacbio)` + `Library prep (Pacbio)` + `Total raw data (Gbp)`, data=filtered_mitVGP)

summary(mod.mitoVGP)
af <- anova(mod.mitoVGP)
af
afss <- af$"Sum Sq"
ml<-cbind(af,PctExp=afss/sum(afss)*100)
ml

write.table(ml, file='../Supplementary Tables/ST4.tsv', quote=FALSE, sep='\t')

mean(mod.mitoVGP$residuals)
acf(mod.mitoVGP$residuals)
lawstat::runs.test(mod.mitoVGP$residuals)

alias(mod.mitoVGP)

ggplot(mod.mitoVGP, aes(`Size selection kbp (Pacbio)`, `Available Pacbio mtDNA reads`)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = `Size selection kbp (Pacbio)`, yend = .fitted), color = "red", size = 0.3)

par(mfrow = c(2, 2))
plot(mod.mitoVGP)

plot(mod.mitoVGP, 4, id.n = 5)

model.diag.metrics <- augment(mod.mitoVGP)
head(model.diag.metrics)

model.diag.metrics <- model.diag.metrics %>%
  mutate(index = 1:nrow(model.diag.metrics)) %>%
  select(index, everything(), -.se.fit, -.sigma)

#model using size selection with 20 kbp cutoff
mod.mitoVGP <- lm(`Available Pacbio mtDNA reads` ~ `Tissue type (Pacbio)` + `Taxonomic group` + `Size selection cutoff` + `DNA extraction (Pacbio)` + `Library prep - fragmentation (Pacbio)` + `Library prep (Pacbio)` + `Total raw data (Gbp)`, data=filtered_mitVGP)

summary(mod.mitoVGP)
af <- anova(mod.mitoVGP)
af
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

