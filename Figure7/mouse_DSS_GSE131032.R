
```{r}
GSE131032_log2_counts_per_million.csv <- read.csv("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE131032_mouse_DSS/GSE131032_log2_counts_per_million.csv.gz", row.names=1)

library(tibble)
```


```{r}
p1 <- GSE131032_log2_counts_per_million.csv[c("Runx1", "Mmp14", "Igf1"),] %>% t() %>% as.data.frame() %>% cbind(data.frame(anno=c(rep(0,3),rep(2,3),rep(4,3),rep(6,3),rep(7,3),rep(8,3),rep(10,2),rep(12,3),rep(14,3)))) %>% rownames_to_column("Sample") %>%  
ggplot(aes(x=anno, y=Runx1)) + 
  geom_point(shape = 21)+
  geom_smooth(alpha = .6, fill="#ffb09c", color="black")+RotatedAxis()+theme(
            axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
        geom_vline(xintercept = c(0,2,4,6,7,8,10,12,14), linetype = 2, color = 'grey')+ 
        guides(color = FALSE, fill = FALSE) + 
        scale_x_continuous(breaks = c(0,2,4,6,7,8,10,12,14))+theme_ArchR()+ggtitle("Runx1")+
     annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = +Inf,
        alpha = .1)



p2 <- GSE131032_log2_counts_per_million.csv["Mmp14",] %>% t() %>% as.data.frame() %>% cbind(data.frame(anno=c(rep(0,3),rep(2,3),rep(4,3),rep(6,3),rep(7,3),rep(8,3),rep(10,2),rep(12,3),rep(14,3)))) %>% rownames_to_column("Sample") %>%  
ggplot(aes(x=anno, y=Mmp14)) + 
  geom_point(shape = 21)+
  geom_smooth(alpha = .6, fill="#ffb09c", color="black")+RotatedAxis()+theme(
            axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
        geom_vline(xintercept = c(0,2,4,6,7,8,10,12,14), linetype = 2, color = 'grey')+ 
        guides(color = FALSE, fill = FALSE) + 
        scale_x_continuous(breaks = c(0,2,4,6,7,8,10,12,14))+theme_ArchR()+ggtitle("Mmp14")+
     annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = +Inf,
        alpha = .1)


p3 <- GSE131032_log2_counts_per_million.csv["Igf1",] %>% t() %>% as.data.frame() %>% cbind(data.frame(anno=c(rep(0,3),rep(2,3),rep(4,3),rep(6,3),rep(7,3),rep(8,3),rep(10,2),rep(12,3),rep(14,3)))) %>% rownames_to_column("Sample") %>%  
ggplot(aes(x=anno, y=Igf1)) + 
  geom_point(shape = 21)+
  geom_smooth(alpha = .6, fill="#ffb09c", color="black")+RotatedAxis()+theme(
            axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
        geom_vline(xintercept = c(0,2,4,6,7,8,10,12,14), linetype = 2, color = 'grey')+ 
        guides(color = FALSE, fill = FALSE) + 
        scale_x_continuous(breaks = c(0,2,4,6,7,8,10,12,14))+theme_ArchR()+ggtitle("Igf1")+
     annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = +Inf,
        alpha = .1)


p4 <- GSE131032_log2_counts_per_million.csv["Cthrc1",] %>% t() %>% as.data.frame() %>% cbind(data.frame(anno=c(rep(0,3),rep(2,3),rep(4,3),rep(6,3),rep(7,3),rep(8,3),rep(10,2),rep(12,3),rep(14,3)))) %>% rownames_to_column("Sample") %>%  
ggplot(aes(x=anno, y=Cthrc1)) + 
  geom_point(shape = 21)+
  geom_smooth(alpha = .6, fill="#ffb09c", color="black")+RotatedAxis()+theme(
            axis.line = element_line(),panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
        geom_vline(xintercept = c(0,2,4,6,7,8,10,12,14), linetype = 2, color = 'grey')+ 
        guides(color = FALSE, fill = FALSE) + 
        scale_x_continuous(breaks = c(0,2,4,6,7,8,10,12,14))+theme_ArchR()+ggtitle("Cthrc1")+
     annotate("rect", xmin = 0, xmax = 7, ymin = -Inf, ymax = +Inf,
        alpha = .1)

cowplot::plot_grid(p1,p2,p3,p4, align = "hv", ncol=5, nrow = 1)




plot_grid(plot1, NULL, plot2, rel_widths = c(1, 0, 1), align = "hv",
          labels = c("A", "B"), nrow = 1)
```



