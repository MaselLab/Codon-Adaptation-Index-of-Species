library(ggplot2)

probfix <- function(x)(1-exp(-0.0005))/(1-exp(-x))

delfix <- function(x)(1-exp(0.0005))/(1-exp(x))

probfixa <- function(x)(1-exp(-0.00005))/(1-exp(-x))

delfixa <- function(x)(1-exp(0.00005))/(1-exp(x))

probfixb <- function(x)(1-exp(-0.000005))/(1-exp(-x))

delfixb <- function(x)(1-exp(0.000005))/(1-exp(x))


plottage <- ggplot(data.frame(x=c(0.000001, 6)), mapping =aes(x)) +
  labs( x = "sN", y = "Ratio of Fixation\ndeleterious vs beneficial\nmutations") +
  theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=28, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))+
  theme(axis.text = element_text(size = 24))+ 
  theme(legend.text = element_text(size = 24, face = "bold"))+ 
  theme(legend.justification=c(1,0.95), legend.position=c(1,0.95))+ 
  theme(legend.title=element_blank())
plottage

plottage <- plottage +stat_function(fun=function(x) delfix(x)/probfix(x), size= 2,linetype = "solid", aes(colour = "s=0.001")) 
plottage <- plottage +stat_function(fun=function(x) delfixa(x)/probfixa(x), size =2,linetype = "dotdash", aes(colour = "s=0.0001")) 
plottage <- plottage +stat_function(fun=function(x) delfixb(x)/probfixb(x), size =2,linetype="dotted",aes(colour = "s=0.00001")) 
plottage + scale_colour_manual(values=c("red", "green", "black")) + scale_linetype_manual(values = c("dotted", "dotdash", "solid"))

plottage1 <- ggplot(data.frame(x=c(0.000001, 6)), mapping =aes(x)) +
  labs( x = "sN", y = "Ratio of Fixation\ndeleterious vs beneficial\nmutations") +
  theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=28, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))+
  theme(axis.text = element_text(size = 24))+ 
  stat_function(fun=function(x) delfix(x)/probfix(x), size= 2,linetype = "solid") 
plottage1
