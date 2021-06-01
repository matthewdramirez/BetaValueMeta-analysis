
# Figure 5 ----------------------------------------------------------------

# This file includes code to reproduce Figure 5 in Ramirez et al. "Meta-analysis of primary producer amino acid δ15N values and their 
# influence on trophic position estimation." 

# Here we use bivariate sensitivy analyses to evaluate how trophic position estimates vary as a function of variation in mean beta values 
# (βGlx-Phe) and trophic discrimination factors (TDFGlx-Phe) for 11 consumers within three model food chains: (A) terrestrial (vascular 
# autotrophs only), (B) freshwater (both vascular and non-vascular autotrophs), and (C) oceanic (non-vascular autotrophs only). 

# Plotted values reflect the difference between a consumer’s TPCSIA estimate for a given β-TDF pairing and the TPCSIA estimate for the baseline 
# β-TDF pairing scenario (stars). The baseline scenario was the TPCSIA estimate derived from the mean TDF for a given panel (7.5 ‰ for primary 
# and secondary consumers, 5.5 ‰ for 3°+ fish consumers, or 4.5 ‰ for 3°+ bird consumers) and a β value of –6.50 (A: terrestrial), –1.75 (B: 
# freshwater), or +3.25 (C: oceanic) that approximate the mean vascular (–6.50) and non-vascular (+3.25) β values resulting from this 
# meta-analysis or their mean (–1.75; i.e., assuming ~ 50% contribution of vascular and non-vascular autotrophs). 

# Isoclines (diagonal black lines) bound bins reflective of a ± 0.25-unit change in trophic position from the baseline scenario TPCSIA estimate. 

# Underlying data were derived from the primary literature and used in conjunction with the single-TDF TPCSIA equation (A, Chikaraishi et al., 
# 2009) or multi-TDF TPCSIA equation (B, C; McMahon & McCarthy, 2016) to estimate consumer TPCSIA. Ranges of β values and TDFs followed known 
# variation for each system and consumer type, such as reductions in mean TDFs with shifts in diet quality or mode of nitrogen excretion (McMahon 
# & McCarthy, 2016).

# Note: We used Adobe Illustrator to compile consumer-specific results panels. 


#### Load Packages

require(plotly)
require(viridis)
#require(stringr)
#require(reshape2)

#> (A) Terrestrial Sensitivity Analysis  ----------------------------------------

### Aphid

# Define beta and TDF values used to estimate TP
Bvec = seq(-12,0,0.25)
TDFvec = seq(6.5,8.5,0.25)

# Create empty matrix to hold TP estimates
aphid.mat = matrix(NA, length(Bvec), length(TDFvec))

# For loop to estimate TP for all combinations of beta and TDF
for(i in 1:length(Bvec)) {   # loop through Betas
  for(j in 1:length(TDFvec)) {    # loop through TDFs
    
    aphid.mat[i,j] <- 1 + ((1.5 - 3.8) - Bvec[i])/TDFvec[j] # Chikaraishi et al. 2009 single-TDF equations
  }
}


aphid.mat
round(max(aphid.mat) - min(aphid.mat),2) # Max-min TP estimate
(TP.base <- aphid.mat[23,5]) # Baseline TP estimate

# Normalize TP estimates to TP baseline
aphid.mat2 <- aphid.mat-TP.base
min(aphid.mat2) #-1.395611
max(aphid.mat2) #1.149843

# Transpose dataset for plotting
t.aphid.mat <- t(aphid.mat2)

# Use plotly to generate heatmap of normalized TP estimates
(aphid <- plot_ly(x= Bvec, y = TDFvec, z = ~t.aphid.mat, type = "contour", 
                  contours = list(showlabels = T, start = -3, end = 3, size = 0.25),
                  colorscale = 'RdBu', reversescale=F,
                  autocontour = F, showscale = F,
                  line = list(color = "black", width=0.5), 
                  width = 425, height = 205)
  %>% layout(title = list(text = "Aphid (TP 2)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

# Export plot as PDF
orca(aphid, file = "aphid.pdf")



### Hoverfly 
Bvec = seq(-12,0,0.25)
TDFvec = seq(6.5,8.5,0.25)

hoverfly.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    hoverfly.mat[i,j] <- 1 + ((9.4 - 3.7) - Bvec[i])/TDFvec[j] # Chikaraishi et al. 2009 single-TDF equations
  }
}


hoverfly.mat
(TP.base <- hoverfly.mat[23,5])
round(max(hoverfly.mat) - min(hoverfly.mat),2)

hoverfly.mat2 <- hoverfly.mat-TP.base
min(hoverfly.mat2) # -0.8385027
max(hoverfly.mat2) #1.213986

t.hoverfly.mat <- t(hoverfly.mat2)


(hoverfly <- plot_ly(x= Bvec, y = TDFvec, z = ~t.hoverfly.mat, type = "contour", 
                     contours = list(showlabels = T, start = -3, end = 3, size = 0.25),
                     colorscale = 'RdBu', reversescale=F,
                     autocontour = F, showscale = F,
                     line = list(color = "black", width=0.5), 
                     width = 425, height = 205)
  %>% layout(title = list(text = "Hoverfly larvae (TP 3)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(hoverfly, file = "hovefly.pdf")


### Parasitoid wasp
Bvec = seq(-12,0,0.25)
TDFvec = seq(6.5,8.5,0.25)

wasp.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    wasp.mat[i,j] <- 1 + ((16.4 - 2.8) - Bvec[i])/TDFvec[j] # Chikaraishi et al. 2009 single-TDF equations
  }
}

wasp.mat
round(max(wasp.mat) - min(wasp.mat),2)
(TP.base <- wasp.mat[23,5])

wasp.mat2 <- wasp.mat-TP.base
min(wasp.mat2) # -1.375318
max(wasp.mat2) #1.847649

t.wasp.mat <- t(wasp.mat2)


(wasp <- plot_ly(x= Bvec, y = TDFvec, z = ~t.wasp.mat, type = "contour", 
                 contours = list(showlabels = T, start = -3, end = 3, size = 0.25),
                 colorscale = 'RdBu', reversescale=F,
                 autocontour = F, showscale = F,
                 line = list(color = "black", width=0.5), 
                 width = 425, height = 205)
  %>% layout(title = list(text = "Parasitoid wasp (TP 4)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(hornet, file = "parasitoid wasp.pdf")


### Hyperparasitoid wasp
Bvec = seq(-12,0,0.25)
TDFvec = seq(6.5,8.5,0.25)

wasp2.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    wasp2.mat[i,j] <- 1 + ((24.0 - 3.2) - Bvec[i])/TDFvec[j] # Chikaraishi et al. 2009 single-TDF equations
  }
}

wasp2.mat
round(max(wasp2.mat) - min(wasp2.mat),2)
(TP.base <- wasp2.mat[23,5])

wasp2.mat2 <- wasp2.mat-TP.base
min(wasp2.mat2) # -1.610526
max(wasp2.mat2) # 2.163636

t.wasp2.mat <- t(wasp2.mat2)


(wasp2 <- plot_ly(x= Bvec, y = TDFvec, z = ~t.wasp2.mat, type = "contour", 
                  contours = list(showlabels = T, start = -3, end = 3, size = 0.25),
                  colorscale = 'RdBu', reversescale=F,
                  autocontour = F, showscale = F,
                  line = list(color = "black", width=0.5), 
                  width = 425, height = 205)
  %>% layout(title = list(text = "Hyperparasitoid wasp (TP 5)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(wasp2, file = "Hyperparasitoid wasp.pdf")

#> (B) Freshwater Sensitivity Analysis -------------------------------------------

require(plotly)
require(viridis)
#require(stringr)
#require(reshape2)

### Caddisfly larvae (TP = 2)

# Define beta and TDF values used to estimate TP
Bvec = seq(-12,6,0.25)
TDFvec = seq(6.5,8.5,0.25)

# Create empty matrix to hold TP estimates
caddisfly.mat = matrix(NA, length(Bvec), length(TDFvec))

# For loop to estimate TP for all combinations of beta and TDF
for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    caddisfly.mat[i,j] <- 1 + ((10.9 - 0.0) - Bvec[i])/TDFvec[j]
    
  }
}

caddisfly.mat
round(max(caddisfly.mat) - min(caddisfly.mat),2)
(TP.base <- caddisfly.mat[42,5])

# Normalize TP estimates to TP baseline
caddisfly.mat2 <- caddisfly.mat-TP.base
min(caddisfly.mat2) #-1.109098
max(caddisfly.mat2) #1.834974

# Transpose dataset for plotting
t.caddisfly.mat <- t(caddisfly.mat2)

# Use plotly to generate heatmap of normalized TP estimates
(caddisfly <- plot_ly(x= Bvec, y = TDFvec, z = ~t.caddisfly.mat, type = "contour", 
                   contours = list(showlabels = TRUE,start = -3, end = 3, size = 0.25),
                   colorscale = 'RdBu', reversescale=F,
                   autocontour = F, showscale = F,
                   line = list(color = "black", width=0.5), 
                   width = 600, height = 205)
  %>% layout(title = list(text = "Caddisfly (TP 2)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

# Export plot as PDF
orca(caddisfly, file = "caddisfly.pdf")



### Trout (TP = 3)

Bvec = seq(-12,6,0.25)
TDFvec = seq(4.5,6.5,0.25)

trout.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
 
        trout.mat[i,j] <- 2 + ((14.7 -  0.0) - 7.5 - Bvec[i])/TDFvec[j] # McMahon and McCarthy 2016 multi-TDF equation
    
  }
}

trout.mat
round(max(trout.mat) - min(trout.mat),2)
(TP.base <- trout.mat[42,5])

trout.mat2 <- trout.mat-TP.base
min(trout.mat2) #-1.052
max(trout.mat2) #1.760308

t.trout.mat <- t(trout.mat2)

(trout <- plot_ly(x= Bvec, y = TDFvec, z = ~t.trout.mat, type = "contour", 
                  contours = list(showlabels = TRUE,start = -3, end = 3, size = 0.25),
                  colorscale = 'RdBu', reversescale=F,
                  autocontour = F, showscale = F,
                  line = list(color = "black", width=0.5), 
                  width = 600, height = 205)
  %>% layout(title = list(text = "Trout (TP 3)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(trout, file = "trout.pdf")


### Osprey

Bvec = seq(-12,6,0.25)
TDFvec = seq(3.5,5.5,0.25)

osprey.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    osprey.mat[i,j] <- 3 + ((20.5 - 5) - 2*7.5 - Bvec[i])/TDFvec[j] # McMahon and McCarthy 2016 multi-TDF equation
    
  }
}

osprey.mat
round(max(osprey.mat) - min(osprey.mat),2)
(TP.base <- osprey.mat[42,5])

osprey.mat2 <- osprey.mat-TP.base
min(osprey.mat2) #-2.068889
max(osprey.mat2) #3.073968

t.osprey.mat <- t(osprey.mat2)


(osprey <- plot_ly(x= Bvec, y = TDFvec, z = ~t.osprey.mat, type = "contour", 
                   contours = list(showlabels = T, start = -3, end = 3, size = 0.25),
                   colorscale = 'RdBu', reversescale=F,
                   autocontour = F, showscale = F,
                   line = list(color = "black", width=0.5), 
                   width = 600, height = 205)
  %>% layout(title = list(text = "Osprey (TP 4)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(osprey, file = "osprey.pdf")





#> (C) Oceanic Sensitivity Analysis --------------------------------------------

require(plotly)
require(viridis)
#require(stringr)
#require(reshape2)

### Copepod (TP = 2)

# Define beta and TDF values used to estimate TP
Bvec = seq(0,6,0.25)
TDFvec = seq(6.5,8.5,0.25)

# Create empty matrix to hold TP estimates
copepod.mat = matrix(NA, length(Bvec), length(TDFvec))

# For loop to estimate TP for all combinations of beta and TDF
for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    copepod.mat[i,j] <- 1 + ((10 - -1.9) - Bvec[i])/TDFvec[j]
    
  }
}

copepod.mat
round(max(copepod.mat) - min(copepod.mat),2)
(TP.base <- copepod.mat[14,5])

# Normalize TP estimates to TP baseline
copepod.mat2 <- copepod.mat-TP.base
min(copepod.mat2) #-0.4917647
max(copepod.mat2) #0.6430769

# Transpose dataset for plotting
t.copepod.mat <- t(copepod.mat2)

# Use plotly to generate heatmap of normalized TP estimates
(copepod1 <- plot_ly(x= Bvec, y = TDFvec, z = ~t.copepod.mat, type = "contour", 
                     contours = list(showlabels = TRUE,start = -3, end = 3, size = 0.25),
                     colorscale = 'RdBu', reversescale=F,
                     autocontour = F, showscale = F,
                     line = list(color = "black", width=0.5), 
                     width = 245, height = 205)
  %>% layout(title = list(text = "Copepod (TP 2)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

# Export plot as PDF
orca(copepod1, file = "copepod.pdf")



### Flying fish (TP = 3)

Bvec = seq(0,6,0.25)
TDFvec = seq(6.5,8.5,0.25)

flyingfish.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    flyingfish.mat[i,j] <- 2 + ((15 - -2.4) - 7.5 - Bvec[i])/TDFvec[j] # McMahon and McCarthy 2016 multi-TDF equation
    
  }
}

flyingfish.mat
round(max(flyingfish.mat) - min(flyingfish.mat),2)
(TP.base <- flyingfish.mat[14,5])

flyingfish.mat2 <- flyingfish.mat-TP.base
min(flyingfish.mat2) #-0.460549
max(flyingfish.mat2) #0.6022564

t.flyingfish.mat <- t(flyingfish.mat2)


(flyingfish <- plot_ly(x= Bvec, y = TDFvec, z = ~t.flyingfish.mat, type = "contour", 
                       contours = list(showlabels = TRUE,start = -3, end = 3, size = 0.25),
                       colorscale = 'RdBu', reversescale=F,
                       autocontour = F, showscale = F,
                       line = list(color = "black", width=0.5), 
                       width = 245, height = 205)
  %>% layout(title = list(text = "Flying fish (TP 3)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(flyingfish, file = "flyingfish.pdf")


### Yellowfin tuna (TP = 4)

Bvec = seq(0,6,0.25)
TDFvec = seq(4.5,6.5,0.25)

yellowfin.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    yellowfin.mat[i,j] <- 2 + ((22.3 - 1.3) - 7.5 - Bvec[i])/TDFvec[j] # McMahon and McCarthy 2016 multi-TDF equation
    
  }
}

yellowfin.mat
round(max(yellowfin.mat) - min(yellowfin.mat),2)
(TP.base <- yellowfin.mat[14,5])

yellowfin.mat2 <- yellowfin.mat-TP.base
min(yellowfin.mat2) #-0.7527273
max(yellowfin.mat2) # 1.087273

t.yellowfin.mat <- t(yellowfin.mat2)


(yellowfin <- plot_ly(x= Bvec, y = TDFvec, z = ~t.yellowfin.mat, type = "contour", 
                      contours = list(showlabels = TRUE,start = -3, end = 3, size = 0.25),
                      colorscale = 'RdBu', reversescale=F,
                      autocontour = F, showscale = F,
                      line = list(color = "black", width=0.5), 
                      width = 245, height = 205)
  %>% layout(title = list(text = "yellowfin tuna (TP 4)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(yellowfin, file = "yellowfin.pdf")



### Swordfish (TP = 5)

Bvec = seq(0,6,0.25)
TDFvec = seq(4.5,6.5,0.25)

swordfish.mat = matrix(NA, length(Bvec), length(TDFvec))

for(i in 1:length(Bvec)) {   #loop through Betas
  for(j in 1:length(TDFvec)) {    #loop through TDFs
    
    swordfish.mat[i,j] <- 2 + ((28.3 - 3.1) - 7.5 - Bvec[i])/TDFvec[j] # McMahon and McCarthy 2016 multi-TDF equation
    
  }
}

swordfish.mat
round(max(swordfish.mat) - min(swordfish.mat),2)
(TP.base <- swordfish.mat[14,5])

swordfish.mat2 <- swordfish.mat-TP.base
min(swordfish.mat2) #-0.8690909
max(swordfish.mat2) #1.255354

t.swordfish.mat <- t(swordfish.mat2)


(swordfish <- plot_ly(x= Bvec, y = TDFvec, z = ~t.swordfish.mat, type = "contour", 
                      contours = list(showlabels = TRUE,start = -3, end = 3, size = 0.25),
                      colorscale = 'RdBu', reversescale=F,
                      autocontour = F, showscale = F,
                      line = list(color = "black", width=0.5), 
                      width = 245, height = 205)
  %>% layout(title = list(text = "Swordfish (TP 5)"))
  %>% layout(xaxis = list(title = "Beta",tickfont = list(size = 20),titlefont = list(size = 20)), 
             yaxis = list(title = "TDF",tickfont = list(size = 20),titlefont = list(size = 20)))
  %>% layout(xaxis = list(dtick = 2, tick0 = 0, tickmode = "linear"), 
             yaxis = list(dtick = 1, tick0 = 0, tickmode = "linear")))

orca(swordfish, file = "swordfish.pdf")

