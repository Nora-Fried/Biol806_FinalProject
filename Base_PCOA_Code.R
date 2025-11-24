library(tidyverse) # loads dplyr, ggplot, and other helpful data processing tools
library(vegan) # for dissimilarity metrics and other ecological functions
library(phyloseq) # Loading just for data
library(ape) # ape has a nice function for comuting pcoas that has a little more functionality than cmdscale() in base R stats package
library(GUniFrac) # if you like unifrac distance rather than bray-curtis you'll need this package
library(ggrepel) # optional - nice to help keep labels on ggplots from running into each other

##installing pyloseq
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

# Load in Global Patterns data from the phyloseq package
data("GlobalPatterns")


#### Step 1: Data cleanup, reformat and rarefy otu table ####
#### ====================================================================== ####
# First let's pull out our metadata for later plotting information
# the conversion to dataframe is necessary because phyloseq objects are full of 
# hidden things that make them incompatible with certain packages. Better to simplify 
# to a dataframe
gp_metadata <- data.frame(sample_data(GlobalPatterns)) %>%
  rename(SampleID = X.SampleID) # sample id column name got messed up so we fix it here

# Extract the taxonomic information - this will be useful for plotting later
gp_tax_info <- data.frame(tax_table(GlobalPatterns))

# Extract OTU table from global patterns dataset and convert it out of phyloseq format
gp_otu_df <- data.frame(otu_table(GlobalPatterns)) 


# Rename the OTUs
# I don't like my OTUs to only be named numbers, so I'll append the word "OTU" in front of them
# this makes it easier so R doesn't mistake them for numbers

# Rename them in the otu table
rownames(gp_otu_df) <- paste0("OTU_", rownames(gp_otu_df))

# Rename them in the taxonomy table
rownames(gp_tax_info) <- paste0("OTU_", rownames(gp_tax_info))

# Also rename them in the phylogenetic tree - not necessary unless you're using something later that requires the phylogenetic tree like unifrac distance does
gp_otu_tree <- phy_tree(GlobalPatterns)
gp_otu_tree$tip.label <- paste0("OTU_", gp_otu_tree$tip.label)

# Transform otu table so species are columns and rows are samples - this is the default that 
# vegan and most other non-microbial ecological packages expects
gp_otu_df_t <- t(gp_otu_df) # note: this converts from dataframe into matrix format


# Figure out rarefaction level
summary(colSums(gp_otu_df)) #min is 58688, this is what we rarefy to below

data.frame(SampleTotCounts = colSums(gp_otu_df)) %>%
  ggplot(aes(x = SampleTotCounts)) +
  geom_histogram()

# This looks like a deeply sampled dataset, it's lowest sample is 58688 K reads,
# Let's try rarefying to that for now.
rar_level = 58688
data.frame(SampleTotCounts = colSums(gp_otu_df)) %>%
  ggplot(aes(x = SampleTotCounts)) +
  geom_histogram() +
  geom_vline(xintercept = rar_level, color = "red")

# Rarefy data:
gp_otu_df_t_rar <- rrarefy(gp_otu_df_t, rar_level)

rowSums(gp_otu_df_t_rar) # check rarefying worked; all should have same read count = rar_level

#### ====================================================================== ####

#### Step 2: Compute distance matrices two ways ####
# different distance metrics say slightly different things about the community,
# it's not necessary to do both, I'm just showing you how for thoroughness 
#### ====================================================================== ####
# Now calculate the distance matrix two ways, bray-curtis distance and unifrac

# Bray-curtis dissimilarities
gp_bc <- vegdist(x = gp_otu_df_t_rar, method = "bray")

# Unifrac dissimilarities 
# note: if your tree is unrooted, unifrac won't work, but midpoint rooting generally is fine
# This might take a while if your tree is long
unifracs <- GUniFrac::GUniFrac(otu.tab = gp_otu_df_t_rar, 
                               tree =  gp_otu_tree, 
                               alpha = c(0, 0.5, 1)) # alpha weighs the importance of lineage abundance; 0 = no weighting, 1 = fully weighted by abundance
gp_uf <- unifracs$unifracs[, , "d_1"] # use the weighted unifrac

#### ====================================================================== ####


#### Step 3: Compute two kinds of ordinations, NMDS and PCoA #####
# generally there's not too much difference between the ordinations, I'm computing them
# here so you know how to do them all
#### ====================================================================== ####
# PCOA oridination, with bray-curtis
gp_bc.pcoa <- pcoa(gp_bc) # this is from the package ape; 
gp_bc.pct_ex <- round((gp_bc.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

# PCOA oridination, with unifrac
gp_uf.pcoa <- pcoa(gp_uf) # this is from the package ape; 
gp_uf.pct_ex <- round((gp_uf.pcoa$values$Relative_eig) * 100, 1) # this line computes the variation on each axis

# Note you can also simply use baseR - see example commented-out example - but the arrow calculations will be a bit trickier
# gp_bc.pcoa <- cmdscale(gp_bc, k = 2, eig = TRUE, add = TRUE) # this line computes the oridnation; K = 2 means only 2 axes worth of points will be shown, you can specify 3 for a 3-d plot, 4 for 4-d, etc
# gp_bc.pct_ex <- round((gp_bc.pcoa$eig/sum(gp_bc.pcoa$eig)) * 100, 1) # this line computes the variation on each axis
# gp_uf.pcoa <- cmdscale(gp_uf, k = 2, eig = TRUE, add = TRUE)
# gp_uf.pct_ex <- round((gp_uf.pcoa$eig/sum(gp_uf.pcoa$eig)) * 100, 1)

# NMDS oridination, with bray-curtis and unifrac
gp_bc.nmds <-  metaMDS(gp_bc, k = 2, trymax = 100) # again, k = # of dimensions - though the impact of choosing k is slightly different in pcoa vs. nmds 
gp_uf.nmds <-  metaMDS(gp_uf, k = 2, trymax = 100) # try upping trymax if no solution is reached


#### ====================================================================== ####

#### Step 5: Prepare plotting dataframes #####
# here we join our sample metadata with the points that were generated by the ordination calculations above
# For simplicity we'll actually just create 1 dataframe with all of the different sets of ordination points
# but you could also make them separately
#### ====================================================================== ####
gp_ord_plot.df <- gp_metadata %>%
  # First join the pcoa data
  left_join(gp_bc.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_BC = Axis.1, PCOA2_BC = Axis.2) %>%
  left_join(gp_uf.pcoa$vectors[,1:2] %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(PCOA1_UF = Axis.1, PCOA2_UF = Axis.2) %>%
  # Next join the nmds data 
  left_join(gp_bc.nmds$points %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(MDS1_BC = MDS1, MDS2_BC = MDS2) %>%
  left_join(gp_uf.nmds$points %>% data.frame() %>%
              rownames_to_column(var = "SampleID"), # here we make a SampleID column before joining the dataframes
            by = "SampleID") %>%
  dplyr::rename(MDS1_UF = MDS1, MDS2_UF = MDS2)


#### ====================================================================== ####

#### Compute Arrows ####
#### ====================================================================== ####
# Learn more about the ins and outs of this here:
# https://www.davidzeleny.net/anadat-r/doku.php/en:suppl_vars_examples

# For PCOA arrows are essentially the correlation between the axis scores and the 
# environmental variable or OTU abundance. For NMDS, it's not quite as simple since
# NMDS's are caluated a little differently

# Note: Arrows assume that there is a linear trend in the variable you are fitting
# across the surface of the plot, this isn't always (or even often) true in unconstrained
# ordination techniques such as PCOA and NMDS, sometimes a better visual representation is 
# to use a more complex visualization of the gradient with ordisurf. 
# See: https://fromthebottomoftheheap.net/2011/06/10/what-is-ordisurf-doing/
# for more information; and https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
# for an example of how to convert ordisurf output into ggplot

# Envfit  - more robust solution works for nmds or PCOA
# this wrapper modifies OTU table to go into envfit, but the same idea can be applied to 
# any numeric variable or a factor that can be converted into a numeric variable. Example: pH, SoilTexture
env_fit_wrapper <- function(ordination, otu_tab, ordination_type = "NMDS", tax_cutoff = NULL) {
  
  # ordination = ordination output
  # ordination_type = either NMDS or PCOA
  # otu_tab = table of organisms abundances in each sample; can be full otu table or subset; columns are organisms
  # tax_cutoff = A number, the abundance cutoff for OTU table, the top xxx most abundant taxa in table will be selected
  
  # Prepare data
  if(ordination_type == "NMDS") {
    ord_axes <- ordination$points
  }
  
  if(ordination_type == "PCOA") {
    ord_axes <- ordination$vectors
  }
  
  # Filter OTU table
  if(!is.null(tax_cutoff)) {
    top_tax <- sort(colSums(otu_tab, na.rm = T), decreasing = T)[1:tax_cutoff]
    names(top_tax)
    otu_filt <- otu_tab[,names(top_tax)]
  } else {
    otu_filt <- otu_tab
  }
  
  # Run envfit
  fit <- envfit(ord = ord_axes, env = otu_filt, perm = 999, na.rm = TRUE)
  
  # adjust p-values for many tests
  pvals.adj <- p.adjust(fit$vectors$pvals, method = "fdr")
  fit$vectors$pval.adj <- pvals.adj
  
  # Use scores to transforms the axes values so the legth of the arrow is propotional to R2
  # Further multiply it for plotting scale with 'OrdiArrowMul"
  scaled_arrow_heads <- data.frame(scores(fit, "vectors"))*ordiArrowMul(fit) 
  
  # Do some data cleaning up for ease of plotting
  arrow_output <- data.frame(scaled_arrow_heads,
                             fit$vectors$arrows,
                             r2 = fit$vectors$r,
                             pvals = fit$vectors$pvals, 
                             pval.adj = fit$vectors$pval.adj) %>%
    rename(Axis.1 = 1, Axis.2 = 2, 
           Axis.1.unscaled = 3, Axis.2.unscaled = 4) %>%
    rownames_to_column(var = "FactorLabel")
  return(arrow_output)
  
} 

gp_bc.nmds.arrows <- env_fit_wrapper(ordination = gp_bc.nmds,
                                     ordination_type = "NMDS",
                                     otu_tab = gp_otu_df_t_rar, 
                                     tax_cutoff = 100) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(gp_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))

gp_uf.nmds.arrows <- env_fit_wrapper(ordination = gp_uf.nmds,
                                     ordination_type = "NMDS",
                                     otu_tab = gp_otu_df_t_rar, 
                                     tax_cutoff = 100)  %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(gp_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))

gp_bc.pcoa.arrows <- env_fit_wrapper(ordination = gp_bc.pcoa,
                                     ordination_type = "PCOA",
                                     otu_tab = gp_otu_df_t_rar, 
                                     tax_cutoff = 100) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but nice
  left_join(gp_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))

gp_uf.pcoa.arrows <- env_fit_wrapper(ordination = gp_uf.pcoa,
                                     ordination_type = "PCOA",
                                     otu_tab = gp_otu_df_t_rar, 
                                     tax_cutoff = 100) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(gp_tax_info %>% 
              rownames_to_column(var = "OTU_ID"),
            by = c("FactorLabel" = "OTU_ID"))

#### ====================================================================== ####

#### Step 5: Plot data #####
#### ====================================================================== ####
# I'm only going to plot a couple of these, but you'll get the idea

# PCOA Bray-Curtis
ggplot(data = gp_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = SampleType)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", gp_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", gp_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = gp_bc.pcoa.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text(data = gp_bc.pcoa.arrows,
            mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel)) +
  # Alternatively
  geom_text_repel(data = gp_bc.pcoa.arrows,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))


# NMDS Bray-Curtis
ggplot(data = gp_ord_plot.df, aes(x = MDS1_BC, y = MDS2_BC)) +
  geom_point(aes(color = SampleType)) +
  # Put scores on axis
  xlab(paste0("MDS1")) + # Can't get % represented in NMDS
  ylab(paste0("MDS2")) + # Can't get % represented in NMDS
  theme_bw() +
  # Arrows
  geom_segment(data = gp_bc.nmds.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text(data = gp_bc.nmds.arrows,
            mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))

#### ====================================================================== ####

# Bonus - Plot by taxonomic group
#### ====================================================================== ####
# Summarize OTU table by class
ClassTab <- t(gp_otu_df_t_rar) %>% 
  data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(gp_tax_info %>% rownames_to_column(var = "OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Class = ifelse(is.na(Class), "ClasssNA", Class)) %>% # fix NA classes
  group_by(Class) %>%
  summarize(across(CL3:Even3, ~sum(.))) %>% # sum class abundances in each sample column
  column_to_rownames(var = "Class") %>%
  t() # retransform data

gp_bc.pcoa.class.arrows <- env_fit_wrapper(ordination = gp_bc.pcoa,
                                           ordination_type = "PCOA",
                                           otu_tab = ClassTab, 
                                           tax_cutoff = 100) %>% 
  # Filter for significant ones only (optional - but helps with plotting)
  filter(pval.adj < 0.05) %>%
  # Bonus, add taxonomy info - not necessary, but  nice
  left_join(gp_tax_info %>% 
              rownames_to_column(var = "OTU_ID") %>%
              select(Phylum:Class) %>% distinct(),
            by = c("FactorLabel" = "Class"))



# Plot Arrows by Class
# PCOA Bray-Curtis
ggplot(data = gp_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = SampleType)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", gp_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", gp_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = gp_bc.pcoa.class.arrows,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = gp_bc.pcoa.class.arrows, max.overlaps = 20,
                  mapping = aes(x = Axis.1, y = Axis.2, label = FactorLabel))
#### ====================================================================== ####

# SuperBonus - Plot an environmental factor instead
#### ====================================================================== ####
# Add a fake numeric and ordinal environmental factor to the metadata
env_data <- gp_metadata %>%
  mutate(Awesomeness = case_when(
    SampleType %in% c("Soil") ~ 10,
    grepl("Sediment", SampleType) ~ 10,
    SampleType %in% c("Feces", "Skin", "Tongue") ~ 3,
    grepl("Freshwater", SampleType) ~ 7,
    grepl("Ocean", SampleType) ~ 5,
    TRUE ~ 1
  ),
  DistanceFromSeaLevel = case_when(
    SampleType %in% c("Soil", "Feces", "Skin", "Tongue") ~ "above",
    grepl("estuary", SampleType) ~ "below",
    grepl("Freshwater", SampleType) ~ "above",
    grepl("Ocean", SampleType) ~ "atsealevel"
  )) %>%
  mutate(DistanceFromSeaLevel = factor(DistanceFromSeaLevel, levels = c("above", "atsealevel", "below"))) %>%  
  select(Awesomeness, DistanceFromSeaLevel)


# Instead of using the wrapper that corrects for p-values we'll use envfit directly
# we only have two columns we're comparing so not correcting our p-values isn't likely to be much of an issue

fit <- envfit(ord = gp_bc.pcoa$vectors, # we're using the Bray-curtis pcoa in this instance 
              env = env_data, perm = 999, na.rm = TRUE)

fit

# Notice how there are both factors and vectors in this fit - that's because we have one 
# numeric column - Awesomeness and another factor column - DistanceFromSeaLevel

# While numeric columns can be represented by arrows, factor levels should be represented by centroids

fit_cont <- (scores(fit, "vectors") * ordiArrowMul(fit)) %>% 
  data.frame() %>%
  rownames_to_column(var = "FactorLabel")

fit_cat <- (scores(fit, "factors") * ordiArrowMul(fit)) %>%
  data.frame() %>%
  rownames_to_column(var = "FactorLabel")


# Plot Arrows by Class
# PCOA Bray-Curtis
ggplot(data = gp_ord_plot.df, aes(x = PCOA1_BC, y = PCOA2_BC)) +
  geom_point(aes(color = SampleType)) +
  # Put scores on axis
  xlab(paste0("Axis 1 (", gp_bc.pct_ex[1],"%)")) + # extract the percent for the first axis
  ylab(paste0("Axis 1 (", gp_bc.pct_ex[2],"%)")) + # extract the percent for the second axis
  theme_bw() +
  # Arrows
  geom_segment(data = fit_cont,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  geom_text_repel(data = fit_cont,
                  mapping = aes(x = Axis.1, y = Axis.2, 
                                label = FactorLabel)) +
  geom_text(data = fit_cat, 
            mapping = aes(x = Axis.1, y = Axis.2, 
                          label = FactorLabel), 
            color = "red")
#### ====================================================================== ####