

```{r}

library(EBImage)

#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fbonly2_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fbonly2_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_noEC2<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))
```




```{r}
#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/Fbonly1_DPAI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/Fbonly1_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_noEC1<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))
```


```{r}

#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fbonly3_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fbonly3_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_noEC3<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))


```

```{r}
#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/25_noTNF_ECs_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/25_noTNF_ECs_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC1<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))
```
```{r}
#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/25fb_ECs_noTNF_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/25fb_ECs_noTNF_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC2<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))
```


```{r}
#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/Fb_ECs3_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/Fb_ECs3_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC3<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))
```

```{r}
#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/Fb_ECs4_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/Fb_ECs4_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC4<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))
```

```{r}

#read in dapi image
dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_1_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_1_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC_TNF1<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))


dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_2_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_2_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC_TNF2<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))


dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_3_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_3_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC_TNF3<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))

dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_4_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_4_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC_TNF4<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))

dna = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_5_DAPI.tiff", type="tiff")
print(dna, short=TRUE)
colorMode(dna) = Grayscale

Runx1 = readImage("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Organoid images/fb_EC_TNF_5_RUNX1.tiff", type="tiff")
print(Runx1, short=TRUE)
colorMode(Runx1) = Grayscale

rgb = rgbImage(blue = dna, red=Runx1)
display(rgb)

nmaskt = thresh(dna, w = 15, h = 15, offset = 0.05)
nmaskf = fillHull( opening(nmaskt, makeBrush(5, shape='disc')) )
dmap = distmap(nmaskf)
nmask = watershed(dmap, tolerance = 10)

display(paintObjects(nmask, rgb, col = "magenta", thick = TRUE))  

st = stackObjects(nmask, rgb)

display(st, all = TRUE)

display( combine(
  toRGB( getFrame(nmaskt, 3) ), 
  colorLabels( getFrame(nmask, 3) )
), all=TRUE )

basic_featu_fb_EC_TNF5<- as.data.frame(computeFeatures.basic(nmask[,,3], Runx1[,,1]))

```


```{r}
```
```{r}
mean(basic_featu_fb_EC1$b.mean)
mean(basic_featu_fb_EC2$b.mean)
mean(basic_featu_fb_noEC1$b.mean)
mean(basic_featu_fb_noEC1$b.mean)

data_fb_EC1<-select(basic_featu_fb_EC1, c(b.mean))
data_fb_EC1$condition<-"EC"

data_fb_EC2<-select(basic_featu_fb_EC2, c(b.mean))
data_fb_EC2$condition<-"EC"

data_fb_EC3<-select(basic_featu_fb_EC3, c(b.mean))
data_fb_EC3$condition<-"EC"

data_fb_EC4<-select(basic_featu_fb_EC4, c(b.mean))
data_fb_EC4$condition<-"EC"


basic_featu_fb_noEC1<-select(basic_featu_fb_noEC1, c(b.mean))
basic_featu_fb_noEC1$condition<-"noEC"

basic_featu_fb_noEC2<-select(basic_featu_fb_noEC2, c(b.mean))
basic_featu_fb_noEC2$condition<-"noEC"

basic_featu_fb_noEC3<-select(basic_featu_fb_noEC3, c(b.mean))
basic_featu_fb_noEC3$condition<-"noEC"

basic_featu_fb_EC_TNF5<-select(basic_featu_fb_EC_TNF5, c(b.mean))
basic_featu_fb_EC_TNF5$condition<-"EC_TNF"
basic_featu_fb_EC_TNF4<-select(basic_featu_fb_EC_TNF4, c(b.mean))
basic_featu_fb_EC_TNF4$condition<-"EC_TNF"
basic_featu_fb_EC_TNF3<-select(basic_featu_fb_EC_TNF3, c(b.mean))
basic_featu_fb_EC_TNF3$condition<-"EC_TNF"
basic_featu_fb_EC_TNF2<-select(basic_featu_fb_EC_TNF2, c(b.mean))
basic_featu_fb_EC_TNF2$condition<-"EC_TNF"
basic_featu_fb_EC_TNF1<-select(basic_featu_fb_EC_TNF1, c(b.mean))
basic_featu_fb_EC_TNF1$condition<-"EC_TNF"



all_data<-rbind(data_fb_EC1, data_fb_EC2, basic_featu_fb_noEC1, basic_featu_fb_noEC2, basic_featu_fb_noEC3, data_fb_EC3, data_fb_EC4)

ggplot(all_data, aes(x=factor(condition, level=c('noEC', 'EC', "EC_TNF")), y=b.mean)) + 
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+ geom_jitter(shape=16, position=position_jitter(0.2))+ theme_classic()+geom_boxplot(width=0.1)



wisteria <- c("grey65", "burlywood3", "khaki2", "plum1", "lightcyan2", "cornflowerblue", "slateblue3")
 
  ggplot(all_data, aes(x=factor(condition, level=c('noEC', 'EC', "EC_TNF")), y=b.mean)) +
  geom_boxplot(aes(fill = condition),            #You can make a box plot too!
               alpha = 0.8, width = 0.7) +      
  geom_point(aes(fill = condition), shape = 21, color = "black", alpha = 0.8,
             position = position_jitter(width = 0.1, seed = 666))+
  scale_fill_manual(values = wisteria[c(3, 1, 6)]) +
  labs(x = "Species",
       y = "Bill length") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(size = 1.2),
        text = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold")
        )




all_data$condition<- as.factor(all_data$condition)


t.test(b.mean ~ condition, all_data, paired=F)
```

```{r}

all_data

summary_data <- all_data %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    Mean = mean(b.mean),
    SD = sd(b.mean)
  )

t_res <- t.test(b.mean ~ condition, data = all_data)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(all_data$b.mean, na.rm = TRUE)


summary_data$condition <- factor(summary_data$condition, levels = c("noEC", "EC"))
all_data$condition <- factor(all_data$condition, levels = c("noEC", "EC"))

summary_data %>% 
ggplot( aes(x = condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, fill=c("red", "lightgrey")) +   geom_jitter(data = all_data, aes(x = condition, y = b.mean), width = 0.15,         # Points
              color = "black", size = 2.5, alpha = 0.6) +            # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size=1) +        # Error bars
    theme_minimal() +
theme(
      panel.background = element_blank(),       # remove inner panel background
      plot.background = element_blank(),        # remove outer background
      panel.grid.major = element_blank(),       # remove major grid lines
      panel.grid.minor = element_blank(),       # remove minor grid lines
      axis.line = element_line(color = "black"), # keep axis lines
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      ) +
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.01,  # space above highest point
           label = sig_label,
           size = 6)
```
