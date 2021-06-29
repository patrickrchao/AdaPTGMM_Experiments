
color_scheme <- data.frame(Method=c(

  "AdaPTGMM",
  "AdaPT","AdaPTg",  "LFDR","Boca Leek",

  "AdaFDR",
  "FDRreg-t","ASH",
  "BH","Storey BH","IHW"),
  col=c(brewer.pal(9,"Set1")[-6],"#000000","#000000",brewer.pal(9,"Set1")[4]),
  lty=c(rep("solid",9),"dashed","dashed"))
color_scheme$lwd <- 2

