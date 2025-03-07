# Default theme set
theme_set(
  theme_bw() +
    theme(
      text = element_text(family = "Arial", size = 12),
      axis.title.x = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold", size = 13),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 10),
    )
)
# Color blind Palette
## 8 colors
cbp8 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
## 12 colors
cbp12 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
## Pair colors
color_yb <- c("#FDB338","#025196") # Yellow/Blue
color_tt <- c("#FFBE6A","#40B0A6") # Tan/Turquoise
color_op <- c("#EB6123","#512888") # Orange/Purple
color_gp <- c("#295E11","#58094F") # Green/Purple
color_br <- c("#2F67B1","#BF2C23") # Blue/Red
color_bp <- c("#10559A","#DB4C77") # Blue/Pink
color_yp <- c("#F4B301","#DB1048") # Yellow/Pink
color_bb <- c("#6A4A3C","#0F65A1") # Brown/Blue
