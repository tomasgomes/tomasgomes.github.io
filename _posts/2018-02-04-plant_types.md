---
title: "What Would Grow in Antarctica?"
layout: post
date: "2018-02-04 21:05:00 +0000"
categories: rblogging
tags: ecology plants random forest
---



Last year an interesting paper [(Wright et al., 2017)](http://science.sciencemag.org/content/357/6354/917) showing how the climate affects leaf size.
Looking attentively at the map in Figure 3, we can see the predictions for leaf size depending on the local climate made by the authors. We see on this map, as well as in those presented in supplementary, that Antarctica is not represented. This likely happens because of lack of data for some parameters and, you know, the lack of plants growing in the continent.

... but what if there were? What could we expect?

To find out, I relied of the dataset published in the [Supplementary Data](http://science.sciencemag.org/content/suppl/2017/08/31/357.6354.917.DC1) section of the paper, where plants are classified in different "Growth Forms" (Trees, Shrubs, Grasses, ...). Let's have a look.

Load libraries for the project

``` r
library(ggplot2) # for plotting
library(ggmap) # for plotting the world map
library(tidyverse) # for some data operations
library(ranger) # for our classifier
library(rgbif) # to retrieve altitude information from Google
## but see https://www.gbif.org/ for more biodiveristy data
```



And now we load the dataset (which I've saved as a csv).
Description for each column can be found in the last sheet of the Excel file.

``` r
leaf_data = read.csv("aal4760-Wright-SM_Data_Set_S1.csv", header = T)
colnames(leaf_data) = gsub(".", "_", colnames(leaf_data), fixed = T)
colnames(leaf_data)
```

    ##  [1] "ID"                     "Site_number"           
    ##  [3] "Site_name"              "Reference"             
    ##  [5] "Genus_species"          "Family"                
    ##  [7] "Order"                  "Name_orig"             
    ##  [9] "TaxonGroup"             "woody_non_woody"       
    ## [11] "Growth_form"            "DecidEver__woody_only_"
    ## [13] "Compound_Simple"        "Leaf_size__cm2_"       
    ## [15] "Whole_leaf_size__cm2_"  "Country"               
    ## [17] "Latitude"               "Longitude"             
    ## [19] "Elevation__m_"          "MAT"                   
    ## [21] "Tgs"                    "TCM"                   
    ## [23] "TCMgs"                  "TWM"                   
    ## [25] "MAP"                    "PPTgs"                 
    ## [27] "cvPPT"                  "MIann"                 
    ## [29] "MIgs"                   "ETq"                   
    ## [31] "ETqgs"                  "RADann"                
    ## [33] "RADgs"                  "RHann"                 
    ## [35] "RHgs"



Let's first see where the measurements were obtained.
Check [this page](https://datadoodlers.blogspot.co.uk/2013/04/r-beginners-plotting-locations-on-to.html) to learn more about plotting maps.

``` r
climate_df = unique(leaf_data[,16:35]) # extract relevant info
world_map <- borders("world", colour="gray50", fill="white") # get map

ggplot()+
  world_map+
  geom_point(data = climate_df, aes(x = Longitude, y = Latitude), 
             colour = "red", size = 1)+
  geom_rug(data = climate_df, aes(x = Longitude, y = Latitude), colour = "black")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "#6faef2"))
```

<img src="plant_types_files/figure-markdown_github/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />



As a side question, we can also look at how do the climate variables relate between these locations

``` r
climate_ann = climate_df[,c(5, 10, 13, 15, 17, 19)] # get only the annual climate variables
climate_ann$dTemp = climate_df$TWM-climate_df$TCM # calculate temperature variation
pca_climate = prcomp(climate_ann, scale. = T, center = T) # PCA of these climate variables

plot(pca_climate) # variance explained by each principal component
```

<img src="plant_types_files/figure-markdown_github/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />


``` r
# prepare data frames for plotting
## PCA points
plot_df = cbind(climate_df[,1:3], pca_climate$x[,1:2])

# PCA vectors for a biplot-like plot
vec_df = data.frame(pca_climate$rotation[,1:2]*8)
vec_df$xst = rep(0, nrow(vec_df))
vec_df$yst = rep(0, nrow(vec_df))
vec_df$name = rownames(vec_df)

ggplot()+
  geom_point(data = plot_df, mapping = aes(x = PC1, y = PC2, colour = Country))+
  geom_segment(data = vec_df, mapping = aes(x = xst, y = yst, 
                                            xend = PC1, yend = PC2), 
               arrow = arrow(angle = 30,length = unit(0.25, "cm"), type = "closed"))+
  geom_text(data = vec_df, aes(x = PC1, y = PC2, label = name), 
            nudge_y = 0.5, fontface = "bold")+
  theme_classic()+
  theme(legend.position = "none",
        aspect.ratio = 1)
```

<img src="plant_types_files/figure-markdown_github/unnamed-chunk-4-2.png" style="display: block; margin: auto;" />



So temperature and light drive PC1, while moisture and percipitation drive PC2. Importantly, we do see some segregation of the sampling locations driven by these variables, meaning that coordinates on a map can be a decent proxy for the climate conditions (At least in this not-so-scientific setting). Unfortunately the legend would be too big and hide the plot, so we won't take this any further.

Now, let's have a look at the Growth Type representation by location (country for simplicity).

``` r
growth_loc = data.frame(table(leaf_data[,c(16, 11)]))

ggplot(growth_loc, aes(x = Country, y = Growth_form, 
                       size = Freq, label = Freq))+
  geom_point()+
  #geom_text()+
  theme_classic()+
  theme(aspect.ratio = 1/2.5, 
        legend.position = "bottom",
        axis.text.x = element_text(size = 5.8, angle = 30, vjust = 1, hjust = 1))
```

![](plant_types_files/figure-markdown_github/unnamed-chunk-5-1.png)

There are some imbalances in the sampling, not necessarily by country (although that could be expected), but some classes are more represented than others.

Let's move on to the predictions! Here. we'll build a model that takes in the Growth Froms as labels and coordinates (Latitude, Longitude and Elevation) as predictors. I also seems fitting that, since we're looking into plants, we would use a Random Forest Classifier (in reality I'm using it because they tend to perform quite well and I'm used to playing with them in single-cell RNA-seq datasets). We'll be using the ranger package for this because it has he fastes implementation of random forests that I know of in R. We will also train two models, a regular one and another where we weight each class according to their abundance.

``` r
# subset the data
sub_leaf = leaf_data[,c(11,17,18,19)]

# regular model
rf_gf = ranger::ranger(Growth_form~., 
                    data = sub_leaf,
                    num.trees = 5000, importance = "impurity", replace = F, 
                    seed = 1)

# adding weights to the matrix (inverse of label frequency)
sub_leaf_w = merge(sub_leaf, 
                   data.frame(1-table(sub_leaf$Growth_form)/nrow(sub_leaf)), 
                   by = 1)
rf_gf_w = ranger::ranger(Growth_form~., 
                    data = sub_leaf_w[,-5],
                    num.trees = 5000, importance = "impurity", replace = F, 
                    seed = 1, case.weights = sub_leaf_w$Freq)

# print importance of the predictors for each model
ranger::importance(rf_gf)
ranger::importance(rf_gf_w)
```

    ##      Latitude     Longitude Elevation__m_ 
    ##     1049.9270      930.4998      587.7785
    ##      Latitude     Longitude Elevation__m_ 
    ##     1101.2313      967.1797      635.4333



Now lets have a look at the confusion matrix for each model. First, the unweighted:

``` r
plot_df = reshape2::melt(apply(rf_gf$confusion.matrix, 1, 
                               function(x) x/colSums(rf_gf$confusion.matrix)))

ggplot(plot_df, aes(x = true, y = predicted, fill = value, 
                    label = round(value, 2)))+
  geom_tile()+
  geom_text(size = 3.2)+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle =30, hjust = 1, vjust = 1))
```

![](plant_types_files/figure-markdown_github/unnamed-chunk-7-1.png)



And now the weighted:

``` r
plot_df2 = reshape2::melt(apply(rf_gf_w$confusion.matrix, 1, 
                               function(x) x/colSums(rf_gf_w$confusion.matrix)))

ggplot(plot_df2, aes(x = true, y = predicted, fill = value, 
                     label = round(value, 2)))+
  geom_tile()+
  geom_text(size = 3.2)+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle =30, hjust = 1, vjust = 1))
```

![](plant_types_files/figure-markdown_github/unnamed-chunk-8-1.png)



The result is (maybe) only marginally better in the weighted model. Interestingly, we see that grasses, herbs and shrubs are among the most confused classes, and these are at least physically similar to each other.

Finally, lets make our preditions. For this I used the elevation function from the rgbif package, which can query the Google Maps API for altitude values based on your cordinates (which I got by expanding a grid). You will need you own access key, which you can find [here](https://developers.google.com/maps/documentation/elevation/start#api_key).

``` r
grid_coord = expand.grid(seq(-90, 90, length.out = 280), 
                         seq(-180, 180, length.out = 280))
colnames(grid_coord) = c("decimalLatitude", "decimalLongitude")

# this query can take a while
map_coord = elevation(input = grid_coord, 
                      key = "AIzaSyDVMGAMHKmyIdpL75h-gxAU44iFxgOKvZw")
colnames(map_coord) = colnames(sub_leaf)[-1]
```



So lets first look at an altitude map. We'll be removing values below 0, since we're not interested about what happens in the ocean.

``` r
map_coord_sub = map_coord[map_coord$Elevation__m_>=0,]

ggplot()+
  world_map+
  geom_point(data = map_coord_sub, size = 0.898,
             aes(x = Longitude, y = Latitude, colour = Elevation__m_))+
  theme_classic()+
  scale_colour_gradientn(colours = terrain.colors(20))+
  theme(panel.background = element_rect(fill = "#6faef2"),
        legend.position = "right")
```

![](plant_types_files/figure-markdown_github/unnamed-chunk-10-1.png)



Looks good! And I was not completely aware of how much Antarctica was above sea level! This can be because of the thinckness of the ice covering the continent (about 1.6km on average), but also because of several moutains and plateaus.

Right, now, on to the preditions using the Random Forests train before.

``` r
# make preditions
pred = predict(rf_gf, data = map_coord_sub)
pred_w = predict(rf_gf_w, data = map_coord_sub)

# add predictions to the map coordinates
map_coord_sub$growth_type = pred$predictions
map_coord_sub$growth_type_w = pred_w$predictions

# format factors to have the same levels
map_coord_sub$growth_type = factor(map_coord_sub$growth_type, 
                                   levels = unique(map_coord_sub$growth_type_w))
map_coord_sub$growth_type_w = factor(map_coord_sub$growth_type_w, 
                                     levels = unique(map_coord_sub$growth_type_w))
```



And we plot the results of the unweighted model:

``` r
ggplot()+
  world_map+
  geom_point(data = map_coord_sub, size = 0.898,
             aes(x = Longitude, y = Latitude, colour = growth_type))+
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  theme_classic()+
  scale_colour_manual(values = RColorBrewer::brewer.pal(8, "Set1"), drop = F)+
  theme(panel.background = element_rect(fill = "#6faef2"),
        legend.position = "bottom")
```

![](plant_types_files/figure-markdown_github/unnamed-chunk-12-1.png)



And those of the weighted one:

``` r
ggplot()+
  world_map+
  geom_point(data = map_coord_sub, size = 0.898,
             aes(x = Longitude, y = Latitude, colour = growth_type_w))+
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  theme_classic()+
  scale_colour_manual(values = RColorBrewer::brewer.pal(8, "Set1"))+
  theme(panel.background = element_rect(fill = "#6faef2"),
        legend.position = "bottom")
```

![](plant_types_files/figure-markdown_github/unnamed-chunk-13-1.png)



They look mostly similar! We can straight up see that these are not great models (unless there are trees buried beneath the sand in the Sahara and crazy lianas in the Australian desert), but they also get some things right (trees in the Amazon, shrubs in southern Argentina and Spain, etc.). I also notice that the region of Mongolia/Northern China has more shrubs and herbs in the weighted model, which should be closer to reality.

Finally, back to our initial question: What would grow in Antarctica? Well, based on this last map, we see that there would be a large divide between herbs in the low lands and shrubs in higher terrain, with trees occupying most of the coastline, allowing for penguins to hide from seals behind them.

Remember that this was an exploratory analysis of this dataset, and followed few biological assumptions. This could probably be seen as a "null model" for plant type distributions, but obviously more rigorous measurements should be used when doing this type of analysis, since terrain coordinates only partially reflect the environment in which plants exist.

Session Info:

``` r
print(sessionInfo(), local = FALSE)
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 16299)
    ## 
    ## Matrix products: default
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] maps_3.2.0      rgbif_0.9.9     ranger_0.9.0    forcats_0.2.0  
    ##  [5] stringr_1.2.0   dplyr_0.7.4     purrr_0.2.4     readr_1.1.1    
    ##  [9] tidyr_0.7.2     tibble_1.3.4    tidyverse_1.2.1 ggmap_2.6.1    
    ## [13] ggplot2_2.2.1  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.13        lubridate_1.7.1     lattice_0.20-35    
    ##  [4] png_0.1-7           assertthat_0.2.0    rprojroot_1.2      
    ##  [7] digest_0.6.12       psych_1.7.8         R6_2.2.2           
    ## [10] cellranger_1.1.0    plyr_1.8.4          backports_1.1.1    
    ## [13] evaluate_0.10.1     httr_1.3.1          RgoogleMaps_1.4.1  
    ## [16] rlang_0.1.4         curl_3.1            lazyeval_0.2.1     
    ## [19] readxl_1.0.0        rstudioapi_0.7      data.table_1.10.4-3
    ## [22] geosphere_1.5-7     whisker_0.3-2       Matrix_1.2-11      
    ## [25] oai_0.2.2           rmarkdown_1.7       urltools_1.7.0     
    ## [28] labeling_0.3        proto_1.0.0         foreign_0.8-69     
    ## [31] triebeard_0.3.0     munsell_0.4.3       broom_0.4.2        
    ## [34] compiler_3.4.2      modelr_0.1.1        pkgconfig_2.0.1    
    ## [37] mnormt_1.5-5        rgeos_0.3-26        htmltools_0.3.6    
    ## [40] geoaxe_0.1.0        crayon_1.3.4        crul_0.5.0         
    ## [43] grid_3.4.2          nlme_3.1-131        jsonlite_1.5       
    ## [46] gtable_0.2.0        magrittr_1.5        scales_0.5.0       
    ## [49] cli_1.0.0           stringi_1.1.5       mapproj_1.2-5      
    ## [52] reshape2_1.4.2      bindrcpp_0.2        sp_1.2-5           
    ## [55] xml2_1.1.1          RColorBrewer_1.1-2  rjson_0.2.15       
    ## [58] tools_3.4.2         glue_1.2.0          hms_0.3            
    ## [61] jpeg_0.1-8          parallel_3.4.2      yaml_2.1.14        
    ## [64] colorspace_1.3-2    rvest_0.3.2         knitr_1.17         
    ## [67] bindr_0.1           haven_1.1.0
