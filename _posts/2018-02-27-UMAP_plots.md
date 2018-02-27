---
title: "Test-driving UMAP on the Tabula Muris dataset"
layout: post
date: "2018-02-27 23:15:00 +0000"
categories: rblogging
tags: [scRNAseq, UMAP, python]
---



Recently I [saw on twitter](https://twitter.com/leland_mcinnes/status/963230617600184320) a reference to a paper in [arXiv](https://arxiv.org/abs/1802.03426) about a new dimensionality reduction method called UMAP (Uniform Manifold Approximation and Projection). Based on the discussion thread and the information presented in the [GitHub page](https://github.com/lmcinnes/umap), it appeared to be a great alternative to tSNE.

tSNE usually gets mixed reactions from people. This is because it is meant to be interpreted as a visualisation tool, even though the clusters it presents are sometimes so well defined that it is tempting to take them (as well as the distances between them) at face value, when in fact it's understood that tSNE cannot be interpreted this way.

On the other hand, UMAP is meant to tackle some of these issues. While I'm not at all equiped to explain how the algorithm works, I understand that it's goal is to learn a manifold that, like tSNE, preserves local distances, but also preserve more of the global structure of the data. Additionally, allowing for the choice of different distance metrics allows for a better embedding of certain data types.

And because tSNE is probably the most used dimensionality reduction method with single-cell RNA-seq data, I decided to try UMAP on such dataset. For this I chose part of the [*Tabula Muris*](https://www.biorxiv.org/content/early/2017/12/20/237446) dataset. Because UMAP is (as far as I know) only implemented in Python, I ran it on a [Jupyter notebook](https://github.com/tomasgomes/UMAP_TM_test). Briefly, after downloading and formatting the data (for which I suggest [this amazing tutorial](https://hemberg-lab.github.io/scRNA.seq.course/tabula-muris.html)), I ran a PCA and passed the 40 most significant PCs to both tSNE and UMAP.

Let's start by loading the libraries we'll need

``` r
library(ggplot2)
library(Laurae) #devtools::install_github("Laurae2/Laurae")
library(cowplot)
library(RColorBrewer)
library(animation)
library(lattice)
library(grid)
library(gridExtra)
```

Now lets load the cell type labels, ad well ad the dimensionality reduction results (see [Jupyter Notebook](https://github.com/tomasgomes/UMAP_TM_test))

``` r
cell_type = read.csv("UMAP/tabula_muris_drop_meta.csv", header = T, row.names = 1)

dimred_list = lapply(list.files("./UMAP/", "mat.csv", full.names = T), 
                     read.csv, header = T, row.names = 1)
names(dimred_list) = substr(list.files("./UMAP/", "mat.csv"), 
                            1, nchar(list.files("./UMAP/", "mat.csv"))-8)
```

One thing that can be seen in the Python notebook is the runtime of tSNE and UMAP. For the same dataset, the sklearn implementation of tSNE took 15-20 minutes to run, whereas UMAPtakes between 1 and 2 minutes at most. UMAP's scalability is made even more evident in the original manuscript.

Now we'll plot all projections, where we have UMAP with a few different distance metrics

``` r
plot_list = list()
for(n in names(dimred_list)){
  plot_df = cbind(dimred_list[[n]][,1:2], cell_type)
  colnames(plot_df) = c("dim1", "dim2", "cell_type", "tissue")
  
  plot_list[[n]] = ggplot(plot_df, aes(x = dim1, y = dim2, colour = tissue))+
    geom_point(size = 0.9, alpha = 0.7)+
    ggtitle(n)+
    guides(colour = guide_legend(ncol = 2))+
    theme_classic()+
    theme(aspect.ratio = 1,
          title = element_text(size = 9),
          legend.text = element_text(size = 7.12),
          legend.key.size = unit(0.35, "cm"))
}
```

Lets look at some of the first plots

``` r
grid_arrange_shared_legend(plot_list$umap_jaccard, plot_list$umap_hamming,
                           plot_list$umap_sokalsneath, plot_list$umap_rogerstanimoto,
                           ncol = 2, nrow = 2,
                           position = "right")
```

<img src="UMAP_plots_files/figure-markdown_github/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

We see that there is no clustering at all with any of these distances. This actually makes sense, since this are all more adequate for binary data.

Let's now look at some more appropriate distances

``` r
grid_arrange_shared_legend(plot_list$umap_canberra, plot_list$umap_correlation,
                           plot_list$umap_cosine, plot_list$umap_minkowski,
                           ncol = 2, nrow = 2,
                           position = "right")
```

<img src="UMAP_plots_files/figure-markdown_github/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

Those look like much clearer clusters! The different metrics do give different resolution and arrangement to the clusters, but overall we can see a nice sepparation of the different mouse tissues, as well as some substructure within these clusters.

Lastly, let's have a look at UMAP run with euclidean distance. Because of the nature of PCA, this distance should be the most adequate to use on those coordinates. We will compare this with the first two principal components and with the tSNE projection

``` r
grid_arrange_shared_legend(plot_list$pca+guides(colour = guide_legend(nrow = 2)), 
                           plot_list$tsne,
                           plot_list$umap_euclidean,
                           ncol = 3, nrow = 1,
                           position = "bottom")
```

<img src="UMAP_plots_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

These look quite nice! And as we can see, UMAP has the clusters with an individual structure mostly similar to tSNE, but in a much clearer disposition. I think that for data **visualisation**, this is clearly a plus.

I just want to highlight some of the subtructure that UMAP gives. Looking at the "Tongue" cluster (I should have probably run it separately, but the goal is also to see how the substructure appears in the global projection)

``` r
plot_df = cbind(dimred_list$umap_euclidean[,1:2], cell_type)
colnames(plot_df) = c("dim1", "dim2", "cell_type", "tissue")
plot_df = plot_df[plot_df$tissue=="Tongue",]

tongue_plot = ggplot(plot_df, aes(x = dim1, y = dim2, colour = cell_type))+
    geom_point(size = 0.9, alpha = 0.7)+
    ggtitle("UMAP Euclidean  - Tongue cell types")+
    guides(colour = guide_legend(ncol = 1))+
    theme_classic()+
    theme(aspect.ratio = 1,
          legend.position = "right",
          title = element_text(size = 9),
          legend.text = element_text(size = 7.12),
          legend.key.size = unit(0.35, "cm"))

print(tongue_plot)
```

<img src="UMAP_plots_files/figure-markdown_github/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

We can clearly see the distinction between clusters of basal cells and keratinocytes. I am not sure about any subclusters of these cell types and I won't comment on these, but this projection is similar to the one found in the original manuscript.

Lastly, we'll have a look at a UMAP 3D projection, just because it's cool (and thanks to [this tutorial](https://rpubs.com/b_t_cooper/lattice_rotation_3D)) on how to make the animated projection

``` r
umap3d = read.csv("./UMAP/umap_euclidean_mat3D.csv", header = T, row.names = 1)

plot_df = cbind(umap3d, cell_type)

cols = brewer.pal(n = 12, name = "Paired")

draw.plot <- function(angles)
{
  for(i in 1:length(angles))
  {
    print(cloud(X0 ~ X1 * X2, data = plot_df, shade = TRUE, distance=0,
        screen=list(z=angles[i],x=-60), groups = tissue,
        par.settings = list(axis.line = list(col = "transparent"),
                            superpose.symbol = list(col = cols, pch = 19)), 
        xlab="",ylab="",zlab=""))
    setTxtProgressBar(txtProgressBar(style=3),i/length(angles))
  }
}

angles = seq(0, 360, 1)
saveGIF(draw.plot(angles), interval = 0.025, ani.width = 550, 
        ani.height = 550, movie.name = "scatter3D.gif")
```

<center>
<img src="UMAP_plots_files/figure-markdown_github/scatter3D.gif" style="display: block; margin: auto;" />
</center>
  
<br />
So, what's the take away?<br />
UMAP appears to provide excelent visualisation capablities, easily separating clusters that usually clump together in a tSNE projection. Improved scalability is also a plus, and the fact that different distance metrics can be used to embbed the data should provide more meaningful data projections. However, UMAP seems to be dependent on several hyperparameters. The major ones are explained in the GitHub tutorial, but the authors have hinted at a future detailed description of their effects (see the twitter thread up top).<br />
