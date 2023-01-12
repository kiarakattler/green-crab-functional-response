library(tidyverse)
theme_paper = function(base_size = 6, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black"),  
      axis.text.x = element_text(size = base_size*1.2, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*1.2, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", size  =  0.2),  
      axis.title.x = element_text(size = base_size*1.5, color = "black", margin = margin(0, 12, 0, 0)),  
      axis.title.y = element_text(size = base_size*1.5, color = "black", angle = 90, margin = margin(0, 5, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "black",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1.2, color = "black"),  
      legend.title = element_text(size = base_size*1.5, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      #panel.border = element_rect(fill = NA, color = "black"),  
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey", color = "grey10"),  
      strip.text.x = element_text(size = base_size*1.5, color = "black"),  
      strip.text.y = element_text(size = base_size*1.5, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

theme_paper_small = function(base_size = 5, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black"),  
      axis.text.x = element_text(size = base_size*1.2, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*1.2, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", size  =  0.2),  
      axis.title.x = element_text(size = base_size*1.5, color = "black", margin = margin(0, 12, 0, 0)),  
      axis.title.y = element_text(size = base_size*1.5, color = "black", angle = 90, margin = margin(0, 5, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "black",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1.2, color = "black"),  
      legend.title = element_text(size = base_size*1.5, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      #panel.border = element_rect(fill = NA, color = "black"),  
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey", color = "grey10"),  
      strip.text.x = element_text(size = base_size*1.5, color = "black"),  
      strip.text.y = element_text(size = base_size*1.5, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

theme_paper_large = function(base_size = 10, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black"),  
      axis.text.x = element_text(size = base_size*1.2, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*1.2, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", size  =  0.2),  
      axis.title.x = element_text(size = base_size*1.5, color = "black", margin = margin(0, 12, 0, 0)),  
      axis.title.y = element_text(size = base_size*1.5, color = "black", angle = 90, margin = margin(0, 5, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "black",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1.2, color = "black"),  
      legend.title = element_text(size = base_size*1.5, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      #panel.border = element_rect(fill = NA, color = "black"),  
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey", color = "grey10"),  
      strip.text.x = element_text(size = base_size*1.5, color = "black"),  
      strip.text.y = element_text(size = base_size*1.5, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

