---
title: "Amateur Myrmecology"
layout: post
date: "2018-06-14 20:09:00 +0000"
categories: rblogging
tags: [ants, SpatialDE, python, R]
---



Myrmecology is the scientific study of ants. I’m not sure I can even be
considered an amateur myrmecologist, but ants are definitely among my
favourite animals.

In my initial blogpost, while looking for a way to obtain elevation
measures given geographical coordinates, I came across this [repository
for biodiversity data](https://www.gbif.org). Naturally, I went to look
for an ant biodiversity dataset, and found the wonderful [Formidabel:
Belgian Ants
Database](https://www.gbif.org/dataset/b528799a-2d52-4023-aa02-9ce081e3ca5f).

In this blog post I will make a very shallow exploration of the dataset,
as well as attempting to find any relevant biological patterns.

## Data Loading and cleanup

Let’s start by loading the libraries needed for this project

``` r
library(ggplot2)
library(ggrepel)
library(ggmap) # for plotting the world map
library(raster)
library(RColorBrewer)
library(animation)
library(cowplot)
library(viridis)
library(rglobi) # get interactions

# a custom theme to plot some maps
map_theme = theme(aspect.ratio = 1,
          title = element_text(size = 8),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank(),
          legend.title = element_blank(), legend.background = element_blank(),
          legend.position = c(0.1,0.25))
```

We will also need to use python later on, so we will setup reticulate to
interface with it.  
Note that, for this to work in RStudio, you currently need
the daily release.

``` r
library(reticulate)
knitr::knit_engines$set(python = reticulate::eng_python)
py_available(initialize = FALSE)
use_python(Sys.which("python"))
py_config()
```

Now lets load the data

``` r
ants_raw = read.table("occurrence.txt", header = T, stringsAsFactors = T, 
                      sep = "\t", comment.char = "", quote = "", row.names = 1)
# ignoring quotes when reading in this dataset is important because of some locality names
# strings as factors help spot variables not to keep
```

And now we proceed to cleaning the data

``` r
ants_filtered = ants_raw[,c("basisOfRecord", "identifiedBy", "samplingProtocol", "eventDate",
                            "habitat", "stateProvince", "locality", "verbatimCoordinates",
                            "decimalLatitude", "decimalLongitude", "genus", "specificEpithet", 
                            "scientificName", "individualCount")]
ants_filtered = ants_filtered[complete.cases(ants_filtered),] # removes samples with NA abundance
ants_filtered = ants_filtered[ants_filtered$samplingProtocol!="",] # removes absence of sampling method
ants_filtered = ants_filtered[ants_filtered$eventDate!="",] # removes absence of sampling date
ants_filtered$eventDate = as.Date(ants_filtered$eventDate)
ants_filtered$habitat[ants_filtered$habitat==""]="Unknown"
```

I also noticed that some entries in the table are actually from the same
capture event. We’re going to add a new variable to be able to track
this, based on the combination of coordinates/date/capture method

``` r
sample_id = apply(ants_filtered[,c("samplingProtocol", "eventDate", "verbatimCoordinates")], 
      1, paste, collapse = "_")
ants_filtered$sample_id = sample_id
```

This leaves us with 5366 sampling events.

## Exploratory Analysis

It will be interesting to look at some variables and see how the data is
distributed according to them. In particular, we should check for
assymetries in the date and method of sample collection, as well as
usage of differen sampling methods in different localities.  
Lets first look at sampling dates:

``` r
plot_df = data.frame(table(ants_filtered$eventDate))
plot_df$Var1 = as.Date(plot_df$Var1)
ggplot(plot_df, aes(x = Var1, y = Freq))+
  geom_point(size = 0.9)+ # the barplot would not show all points
  scale_x_date(date_breaks = "10 years", date_labels = "%Y")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust =1),
        aspect.ratio = 1/1.75)
```

<img src="ants_files/figure-markdown_github/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

Most samples were collected since about 1997, so we will be filtering by
records starting in 1997-01-01 to have a sense of the corrent ant
distribution. As a sidenote, although not depicted, removing records
with non-registered capture method removed most of the older records.  
Now it’s time to look into sampling methods’ frequency and capture
efficiency

``` r
# filtering by data
ants_filtered = ants_filtered[ants_filtered$eventDate>="1997-01-01",]

# looking into sampling methods
plot_df = data.frame(table(ants_filtered$samplingProtocol))
plot_df$Var1 = as.character(plot_df$Var1)
sub_data = ants_filtered[,c("sample_id", "samplingProtocol", "individualCount")]
sub_data = aggregate(individualCount ~ sample_id+samplingProtocol, data = sub_data, FUN=sum)
sub_data$samplingProtocol = as.character(sub_data$samplingProtocol)
sub_data1 = aggregate(individualCount ~ samplingProtocol, data = sub_data, FUN=sum)
sub_data2 = aggregate(individualCount ~ samplingProtocol, data = sub_data, FUN=mean)
plot_df = Reduce(function(x, y) merge(x, y, by = 1), 
                 list(plot_df, sub_data1, sub_data2), accumulate=FALSE)
colnames(plot_df) = c("samplingProtocol", "Freq", "Counts", "Mean")

ggplot(plot_df, aes(x = Counts, y = Freq, 
                    label = samplingProtocol, size = Mean))+
  geom_point()+
  geom_label_repel(size = 2.85, label.padding = unit(0.1, "cm"),
                   min.segment.length = 0)+
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000))+
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000))+
  theme_classic()+
  theme(aspect.ratio = 1/1.15)
```

<img src="ants_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

This shows a clear preference for Pitfall and Manual sampling methods.
Additionally, it reveals that Yellow color traps and White pitfalls have
high capture rates. Indeed, [it has been
shown](https://www.eje.cz/artkey/eje-201002-0018_effect_of_the_colour_of_pitfall_traps_on_their_capture_efficiency_of_carabid_beetles_coleoptera_carabidae_s.php)
that ants (as well as other arthropods) have a preference for these
colours, however a further breakdown by environment would be needed to
reach a better conclusion.

Finally, before looking at the association between ant species, lets
observe how sampling has changed in the last 20 years (see [this
page](https://stackoverflow.com/questions/26914616/r-how-to-add-border-countries-to-a-country-spatialpolygons-map)
for reference)

``` r
# download country map references (raster package)
belgium <- getData("GADM",country="Belgium",level=2)
ned <- getData("GADM",country="Netherlands",level=0)
france <- getData("GADM",country="France",level=0)
ger <- getData("GADM",country="Germany",level=0)
lux <- getData("GADM",country="Luxembourg",level=0)


map_bel = ggplot(belgium,aes(x=long,y=lat,group=group))+
  geom_polygon(color="grey37", fill="grey82")+
  geom_polygon(data=ned, fill="grey58",color="grey41")+
  geom_polygon(data=france, fill="grey58",color="grey41")+
  geom_polygon(data=ger, fill="grey58",color="grey41")+
  geom_polygon(data=lux, fill="grey58",color="grey41")+
  coord_map(xlim=c(-0.1,.1)+bbox(belgium)["x",],ylim=c(-0.1,.1)+bbox(belgium)["y",])+ #centre image
  theme_classic()+
  theme(panel.background = element_rect(fill = "#6faef2"))

draw.maps = function(df){
  # colours for habitats
  pal = c("#bfed42", "#f7d200", "#26efca", "#1b8732", "#929392",
          "#db4f48", "#edad6d", "#61ce66", "#77550b", "#538aaa")
  names(pal) = unique(df$habitat)
  df$habitat = factor(df$habitat, levels = unique(df$habitat))
  
  df$bins = cut(df$eventDate, "6 months")# make 6mo intervals
  
  # date dataframe
  date_df = data.frame(p = 1:length(levels(df$bins)), date = levels(df$bins))
  
  for(d in levels(df$bins)){
    # date plot
    date_df$high = ifelse(date_df$date==d, "bold", "plain")

    date_plt = ggplot(date_df, aes(x = p, y = "Date", label = date))+
      geom_text(angle = 40,  
                mapping = aes(fontface = high, colour = high, size = high))+
      scale_size_discrete(range = c(4.15, 3.25))+
      scale_colour_manual(values = c("black", "grey65"))+
      labs(y = "Date")+
      scale_x_continuous(breaks = seq(1, 32, 4))+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0, 0.25, 0), "cm"),
            aspect.ratio = 1/6)
    
    # map
    sub_df = df[df$bins==d,]
    sub_df = aggregate(individualCount ~ habitat+decimalLatitude+decimalLongitude, 
                       data = sub_df, FUN=sum)
    
    plt_samp = map_bel+
      geom_point(data = sub_df, mapping = aes(x = decimalLongitude, y = decimalLatitude,
                                              group = habitat,
                                              colour = habitat, size = log2(individualCount+1)))+
      scale_colour_manual(values = pal, drop = F)+
      scale_size_continuous(limits = c(1, 12))+
      guides(colour = guide_legend(order = 1), size = guide_legend(order = 2))

    joint_plot = plot_grid(date_plt, plt_samp, nrow = 2, rel_heights = c(0.35,1))
    print(joint_plot)
  }
  sub_df = aggregate(individualCount ~ habitat+decimalLatitude+decimalLongitude, 
                     data = df, FUN=sum)
  
  plt_samp = map_bel+
    geom_point(data = sub_df, mapping = aes(x = decimalLongitude, y = decimalLatitude,
                                            group = habitat,
                                            colour = habitat, size = log2(individualCount+1)))+
    scale_colour_manual(values = pal, drop = F)+
    scale_size_continuous(limits = c(1, 12))+
    guides(colour = guide_legend(order = 1), size = guide_legend(order = 2))+
    ggtitle("All Sites")
  print(plt_samp)
  print(plt_samp)
  print(plt_samp)
  
}
saveGIF(draw.maps(ants_filtered), interval = 1.25, ani.width = 550,
        ani.height = 550, movie.name = "samplesTime.gif")
```

<center>
<img src="ants_files/figure-markdown_github/samplesTime.gif" style="display: block; margin: auto;" />
</center>
<br />
We also notice in this map that there is some bias towards sampling in
the Flanders region. To avoid any bias created by the “Wallonian
outliers”, provinces from this region will be removed from further
analysis.

``` r
ants_filtered = ants_filtered[ants_filtered$stateProvince %in% c("Limburg", "Oost-Vlaanderen",
                                                                 "Vlaams-Brabant", "Antwerpen",
                                                                 "West-Vlaanderen",
                                                                 "Brussel Hoofdstedelijk G"),]

# and removing samples without names
ants_filtered = ants_filtered[ants_filtered$genus!="",]
```

## Ant species spatial association

Different ant species can compete (or not) for the same resources. We
can then ask which ant species follow similar distribution patterns.  
You should have realised by this point that I am a (very) amateur
ecologist, and am for the first time playing with this type of analysis.
Intuitively, I will start by using something familiar. My friend and
former lab mate [Valentine](https://twitter.com/vallens) has recently
developed a [method](https://github.com/Teichlab/SpatialDE) for spatial
transcriptomics capable of detecting spatial patterns of gene
expression. But how will it work for ants?  
Lets start by formatting our data

``` r
metadata = ants_filtered[,c("decimalLatitude", "decimalLongitude", "individualCount")]
metadata = aggregate(individualCount~decimalLatitude+decimalLongitude, 
                     data = metadata, FUN = sum)
rownames(metadata) = paste0("c_",metadata$decimalLatitude, "x", metadata$decimalLongitude)

species_data = ants_filtered[,c("decimalLatitude", "decimalLongitude", 
                                "individualCount", "genus", "specificEpithet")]
species_data$spec = paste(species_data$genus, species_data$specificEpithet)
species_data$sample = paste0("c_",species_data$decimalLatitude, "x", species_data$decimalLongitude)
species_data = species_data[,c(6,7,3)]
species_data = reshape2::acast(species_data, formula = spec~sample, 
                               value.var = "individualCount", fun.aggregate = sum)
species_data = data.frame(species_data[,rownames(metadata)])

# filter out low counts
species_data = species_data[rowSums(species_data>0)>2,]

# plot sample locations with total number of individuals captured
map_bel+
  geom_point(data = metadata, mapping = aes(x = decimalLongitude, y = decimalLatitude, 
                                            colour = log2(individualCount)), group = 1)+
  ggtitle("log2(counts) per site")+
  map_theme
```

<img src="ants_files/figure-markdown_github/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

And now we can go to python.

``` python
# Load packages
import SpatialDE
import NaiveDE
import pandas as pd
import numpy as np
np.random.seed(1)
# prepare data
## Anscombe normalization (normal distributed data)
norm_expr = NaiveDE.stabilize(r.species_data).T
## regress out total counts
resid_expr = NaiveDE.regress_out(r.metadata, norm_expr.T, 'np.log2(individualCount)').T
# Run SpatialDE
X = r.metadata[['decimalLatitude', 'decimalLongitude']]
results = SpatialDE.run(X, resid_expr)
```

Now lets plot 2 of the top and bottom ant species

``` r
yyy = merge(metadata, py$norm_expr, by = 0)
rownames(yyy) = yyy[,1]
yyy = yyy[,-1]
colnames(yyy) = gsub(pattern = " ", replacement = "_", colnames(yyy))

test_df = py$results[py$results$l>=0.1,]
test_df = test_df[order(test_df$qval, test_df$l, decreasing = c(F,T)),]

l_plots = list()
for(n in test_df$g[c(1,2,(nrow(test_df)-1), nrow(test_df))]){
  vals = paste0(n, " (l = ", round(test_df[test_df$g==n,"l"], 2),
               ", qval = ", round(test_df[test_df$g==n,"qval"], 2),")")
  n = gsub(" ", "_", n)
  l_plots[[n]] = map_bel+
  geom_point(data = yyy, mapping = aes_string(x = "decimalLongitude", y = "decimalLatitude",
                                                  colour = n), group = 1, size = 1)+
    scale_color_viridis(option = "magma", limits = c(-4,7))+
    ggtitle(vals)+
    guides(colour = guide_colourbar(barwidth = 1))+
    map_theme
}
cowplot::plot_grid(plotlist = l_plots)
```

<img src="ants_files/figure-markdown_github/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

We can clearly see that the top two ants are present in specific areas, whereas the
bottom two are either completely dispersed or present in very low numbers.

Now based on the ants with spatial variation, we should be able to find
those with similar distributions using… uh… “Automated Expression
Histology”?

``` python
sign_results = results.query('qval <= 0.05')
# we'll choose 3 clusters and use approximately the mean lengthscale of the significant results
np.random.seed(1)
histology_results, patterns = SpatialDE.aeh.spatial_patterns(X, resid_expr, sign_results, 
C=3, l=0.22, verbosity=1)
patterns.columns = ["cl0","cl1","cl2"]
```

And now we can plot the realizations of the ant distribution clusters

``` r
cldf = merge(py$patterns, metadata, by = 0)

l_plots = list()
for(n in paste0("cl", 0:2)){
  l_plots[[n]] = map_bel+
  geom_point(data = cldf, mapping = aes_string(x = "decimalLongitude", y = "decimalLatitude",
                                              colour = n), group = 1, size = 1)+
    scale_color_viridis(option = "magma", limits = c(-0.8,1.8))+
    ggtitle(n)+
    guides(colour = guide_colourbar(barwidth = 1))+
    map_theme
}
cowplot::plot_grid(plotlist = l_plots, ncol = 3)
```

<img src="ants_files/figure-markdown_github/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

Finally, we can look for annotated interactions between the clustered
species in [GLOBI](https://www.globalbioticinteractions.org/). The
spatial clustering is able to recapitulate some previously known
interactions. This could also be compared with information about ant
behaviour to establish coexistance relationships between different
species.

``` r
ant_int_list = list()
for(ant in py$histology_results$g){
  ant_int_list[[ant]] = unique(get_interactions(taxon = ant, 
                                                interaction.type = "interactsWith")[,c(2,5,7)])
  ant_int_list[[ant]] = ant_int_list[[ant]][ant_int_list[[ant]][,3] %in% py$histology_results$g,]
}

int_mat = unique(Reduce(rbind, ant_int_list)[,c(1,3)])
int_mat = reshape2::dcast(int_mat, source_taxon_name ~ target_taxon_name, drop = F)
rownames(int_mat) = int_mat[,1]
int_mat = int_mat[,-1]
for(i in  1:nrow(int_mat)){
  for(j in 1:ncol(int_mat)){
    int_mat[i,j] = if(is.na(int_mat[i,j])) "0" else "1"
  }
}
int_mat = apply(int_mat, 1, as.numeric)
rownames(int_mat) = colnames(int_mat)
meta = data.frame(row.names = py$histology_results[,1],
                  "cluster" = paste0("cl",py$histology_results[,3]))

pheatmap::pheatmap(int_mat, annotation_col = meta, annotation_row = meta, 
                   color = c("grey70", "red"), legend = F,
                   clustering_distance_rows = "manhattan",
                   clustering_distance_cols = "manhattan",
                   treeheight_col = 20, treeheight_row = 20)
```

<img src="ants_files/figure-markdown_github/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

One of the largest drawbacks of this analysis is the geographical range 
of our samples. While we might be able to associate some species and
interactions with broad regions (for instance, ants in the coastal region),
the habitat diversity at a subcluster level can stop us from establishing
deeper associations between ant communities and their environment.

