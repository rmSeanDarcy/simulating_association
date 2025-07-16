###########################################################################################
### SAMSARA - Loading Set of analyses data and preparing for plotting                   ###
###########################################################################################
#options(bitmapType = "Xlib")

###########################################################################################
### ggplot2 plotting modifications
# Control font sizes
sml_font <- 8
med_font <- 10
big_font <- 12
font_size_control <- theme(text=element_text(size=sml_font), #change font size of all text
                           axis.text=element_text(size=sml_font), #change font size of axis text
                           axis.title=element_text(size=med_font), #change font size of axis titles
                           plot.title=element_text(size=big_font), #change font size of plot title
                           legend.text=element_text(size=med_font), #change font size of legend text
                           legend.title=element_text(size=med_font)) #change font size of legend title   

###########################################################################################
### Save as .svg wrapper:
save_svg <- function(file, fig, width, height, pointsize = 10) {
  # Convert cm to inches for svglite
  width_in <- width / 2.54
  height_in <- height / 2.54
  svglite::svglite(file, width = width_in, height = height_in, pointsize = pointsize)
  print(fig)
  dev.off()
}
# save_svg <- function(file, fig, width, height, pointsize = 10) {
#   # Convert cm to inches for svglite
#   width_in <- width / 2.54
#   height_in <- height / 2.54
#   Cairo::CairoSVG(filename = file, width = width_in, height = height_in, pointsize = pointsize)
#   print(fig)
#   dev.off()
# }


# Remove axes and legend if necessary
cxy <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),
             axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())
cxyl <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),
              axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")
cx <- theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())
cxl <- theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.position = "none")
cy <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())
cyl <- theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank(),legend.position = "none")
cl <- theme(legend.position = "none")

###########################################################################################
### Get info about quartiles in BC resource preference data
get_Qs <- function(bcdm, r_num) {
  mean_Q1 <- c()
  mean_Q3 <- c()
  for (r in r_num) {
    bcdmr <- bcdm[bcdm$repl == r,] 
    mean_Q1 <- c(mean_Q1, quantile(as.numeric(bcdmr[lower.tri(bcdmr)]), c(.25)))
    mean_Q3 <- c(mean_Q3, quantile(as.numeric(bcdmr[lower.tri(bcdmr)]), c(.75)))
  }
  bc_thr_neg <- mean(mean_Q1)           # Lower threshold for testing predictions
  bc_thr_pos <- mean(mean_Q3)           # Upper threshold for testing predictions
  res <- c(bc_thr_neg,bc_thr_pos)
  return(res)
}

###########################################################################################
### Plotting colours and parameters:
int_col <- '#ffa500ff'    #'#F1DAFFFF'
eucint_col <- '#008080ff'
euc_col <- '#008080ff'
rest_col <- '#D1D3D2FF'
multipl_q = (((10*10)-10)/2)*0.25
multipl_posint = 4
multipl_negint = 16
boxalpha = 0.8
boxcolor = 'white'
titltxt <- 'blank'
hnum <- 25


