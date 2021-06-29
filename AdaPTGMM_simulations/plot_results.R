library(gridExtra)
library(grid)

library(tools)
plot_results <- function(full_log,testing){


  testing <- gsub("_"," ",testing)
  filename_title <- paste(toTitleCase(testing),"Testing")

  col <- as.character(color_scheme$col)
  names(col) <- as.character(color_scheme$Method)
  lty <- as.character(color_scheme$lty)
  names(lty) <- as.character(color_scheme$Method)
  line_size <- 1

  old_log <- full_log
  full_log <-
    full_log %>% group_by(Method, Alpha) %>%  summarise_at(.vars=c("TPR","FDR","Rej"), funs(mean(., na.rm=TRUE)))

  file_name <- paste(filename_title,"FDR")
  plot_title <- "False Discovery Rate"

  full_log$FDR_diff <- full_log$FDR - full_log$Alpha
  p1 <-
    full_log %>% filter(Alpha < 0.301 & Alpha > 0) %>%
    ggplot(aes( x = Alpha,
                y = FDR_diff,
                group = Method,
                color = Method,
                linetype = Method,
                linewidth = Method,
                fill = Method
    )) + annotate("rect",xmin=-Inf,xmax = Inf, ymin = 0, ymax = Inf,alpha=0.3,fill="pink")+
    geom_abline(
      intercept = 0,
      slope = 0,
      alpha = 0.7,
      linetype="solid",
      color = "gray40",
      size=0.5
    ) +
    geom_line(show.legend = TRUE, alpha = 0.8,size=line_size) +
    labs(title = plot_title)  +
    xlab(expression(paste("Target FDR Level ", alpha))) + ylab(expression(paste("FDR " -  alpha,"  Level"))) + theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 25),
      legend.position = "none"
    ) +
    guides(fill = guide_legend(keywidth = 0.15, keyheight = 0.15)) + scale_linetype_manual(values = lty) +
    scale_color_manual(values = col)+ scale_x_continuous(limits=c(0.02,0.2),expand = c(0.01, 0.01)) +
    scale_y_continuous(limits=c(-0.2,0.05), expand = c(0, 0))



  ggsave( paste0("Images/", file_name, ".png"), width = 20, height = 15, dpi = 200, units = "cm")



  file_name <- paste(filename_title,"Rejections")
  plot_title <- "Number of Rejections"

  p2 <-
    full_log %>% filter(Alpha < 0.301 & Alpha > 0) %>% ggplot(aes(
      x = Alpha,
      y = Rej,
      fill = Method,
      color = Method,
      linetype = Method,
      linewidth = Method

    )) +
    geom_line(show.legend = TRUE, alpha = 0.8,size=line_size) +
    labs(title = plot_title) + theme_classic() +
    xlab(expression(paste("Target FDR Level ", alpha))) + ylab("Number of Rejections") +
    theme(
      text = element_text(size = 20),
      legend.position = "none",
      plot.title = element_text(size = 25)
    ) + scale_y_continuous( limits=c(0,NA),expand = c(0, 0))+
    scale_linetype_manual(values = lty) + scale_color_manual(values = col)+ scale_x_continuous(expand = c(0.01, 0.01))


  ggsave( paste0("Images/", file_name, ".png"), width = 20, height = 15, dpi = 200, units = "cm")


  file_name <- paste(filename_title,"TPR")
  plot_title <- "True Positive Rate"
  p3 <-
    full_log %>% filter(Alpha < 0.301 & Alpha > 0) %>% ggplot(aes(
      x = Alpha,
      y = TPR,
      fill = Method,
      color = Method,
      linetype = Method,
      linewidth = Method,
    )) +
    geom_line(show.legend = TRUE, alpha = 0.8,size=line_size) +
    labs(title = plot_title) + theme_classic() +
    xlab(expression(paste("Target FDR Level ", alpha))) + ylab("TPR") +
    theme(
      text = element_text(size = 20),
      legend.background = element_rect(fill = 'transparent'),
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      plot.title =  element_text(size = 25),
      legend.key.size=unit(1,"cm")
    ) + scale_y_continuous( limits=c(0,NA),expand = c(0, 0)) +
    scale_linetype_manual(values = lty) + scale_color_manual(values = col)+ scale_x_continuous(expand = c(0.01, 0.01))
  #print(p3)

  ggsave(paste0("Images/", file_name, ".png"), width = 20, height = 15, dpi = 200,  units = "cm")
  #
  #   print(full_log %>% filter(Alpha<0.101 & Alpha > 0)%>%ggplot(aes(x=Alpha,y=TPR,fill=Method,color=Method))+geom_line(show.legend=TRUE,alpha=0.8) +
  #           labs(title = plot_title)+theme(text = element_text(size=20))+ylim(0,1)+
  #           xlab(expression(paste("Target FDR Level ",alpha)))+ylab("TPR"))
  #
  #   ggsave(paste0("Images/",file_name,"_Zoom",".png"),width=20,height=15,dpi=200,units="cm")




  #full <- arrangeGrob(p1,p3,ncol=2)
  #print(full)

  print("@@@@@")
  #print(out)
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  out <- grid.arrange(arrangeGrob(p1, p3+theme(legend.position="none"),nrow=1),g_legend(p3), nrow = 1,widths=c(6,1.5),
                      top = textGrob(filename_title,gp=gpar(fontsize=30)),padding=unit(1, "line"))
  file_name <- paste(filename_title,"full_plots",sep="_")
  ggsave(paste0("Images/", file_name, ".pdf"), out, width = 40, height = 15, dpi = 200, units = "cm")


}
# shiny_plot_results(readRDS("data/final/High_Dimensional One Sided.RDS"),"Logistic One Sided")
# shiny_plot_results(readRDS("data/final/High_Dimensional Two Sided.RDS"),"Logistic Two Sided")
# shiny_plot_results(readRDS("data/final/High_Dimensional Interval.RDS"),"Logistic Interval")
#
#shiny_plot_results(readRDS("data/combined/Logistic One Sided.RDS"),"Logistic One Sided")
#shiny_plot_results(readRDS("data/combined/Logistic Two Sided.RDS"),"Logistic Two Sided")
#shiny_plot_results(readRDS("data/combined/Logistic Interval.RDS"),"Logistic Interval")


#shiny_plot_results(readRDS("data/combined_exp4_full/Logistic One Sided.RDS"),"Logistic One Sided")
#shiny_plot_results(readRDS("data/combined_exp4_full/Logistic Two Sided.RDS"),"Logistic Two Sided")
#shiny_plot_results(readRDS("data/combined_exp4_full/Logistic Interval.RDS"),"Logistic Interval")
