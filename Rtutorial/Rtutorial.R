##########################################
# Name: MacIntosh Cornwell
# Email: macintosh.cornwell@nyulangone.org
##########################################
## This is an R Tutorial used for the purposes of demonstrating basic R functions, in addition to a few modules that work with data manipulation, plotting, and differential expression

## Load in Libraries
#packagelist = c()
#junk = lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

## One of the first things we want to do is read in our data to play with - or obtain it from some sort of online resource
## Inquiring via GEO to obtain the supplementary file for this experiment: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94493
# library("GEOquery")
# getGEOSuppFiles(GEO = "GSE94493", makeDirectory=TRUE, baseDir = "/Users/tosh/Desktop/Ruggles_Lab/code/GEO_output/", fetch_files = TRUE)

############
# START HERE
############
##### BASIC READING IN OF DATA AND INDEXING OUR DATA
## This Rscript is designed to be run through piece by piece - and to later use as a reference for further analyses. single lines of this script can be run in the shell by pressinc COMMAND+ENTER (mac) or CONTROL+ENTER(pc)

## Read in the local files
# Define the file that we are reading in:
intablefile = "/Users/mgc439/Code/codingclub/Rtutorial/GSE94493_ER_Mutants_Cufflinks_Gene_Counts.csv"
# Read in the file as a table, our file has a HEADER - so set this to true, and we want to use the FIRST column as ROW NAMES - so we set row.names=1, and the separator for this file is a COMMA (its a CSV - comma separated values) (NOTE - if this was a tab separated file (txt) then the separator is a "\t")
intable = read.table(intablefile, header=TRUE, row.names = 1, sep=",")


# Repeat for the metadata - NOTE though that we arent assigning row names here
# we are also adding a flag here - stringsAsFactors = FALSE - this is a little outside the scope of a tutorial, but for now, make sure to include this when reading in tables so that all of your strings as read as strings!! (maybe more on this later..)
inmetafile = "/Users/mgc439/Code/codingclub/Rtutorial/rtutorial_metadata.csv"
metatable = read.table(inmetafile, header=TRUE, sep=",", stringsAsFactors = FALSE)


## What are we looking at?
# dim(<object>) - returns the dimensions of the table
dim(intable) # 20111, 69
dim(metatable) # 48,5
#head(<object>) - returns the top 6 lines (and all columns!!!) of the table
head(intable)
head(metatable)


## Out data table has way too many columns (69) and we dont need to see all of them right now, we just want a snapshot of the data - so lets look at just the first 5 rows, and the first 10 columns:
# <TABLE>[ROWS,COLUMNS]
intable[1:5,1:10] # rows 1-5 and columns 1-10
metatable[1:5,] # rows 1-5 and ALL COLUMNS


## Now that we know how to index the data using numbers - lets try something else.. indexing by strings!
## So first of all - run the dimensions of your data again, note how there are 69 columns of data, but 72 rows of metadata
# Now we want to subset our metadata to only include the samples that we have in our datatable - there are a few ways we can do this, but lets start with basics
# First - you can index your data by row names and colnames as well as numbers
intable["FAM173A", c("IRJ1", "IRJ2")] # FAM173A row for just samples IRJ1 and IRJ2
metatable[1,"sample"] # IRJ1


# Next - the row and colnames for objects can be called as objects unto themselves
head(colnames(intable)) # "IRJ1" "IRJ2" "IRJ3" "IRJ5" "IRJ6" "IRJ7"


# So - if colnames(intable) returns all of the column names - which are the sample IDs for this dataset - how can you use this to only see which entries of metadata that we have that are NOT in our intable
# so we want to know which of our COLNAMES of TABLE are IN our FIRST COLUMN of METATABLE
colnames(intable) %in% metatable[,1] # Returns vector of 69 entries - some true and some false
# what about the opposite
metatable[,1] %in% colnames(intable) # returns a vector of 48 entries - 2 of which are FALSE
#NOTE that TRUE/FALSE can also be intpretted as 1s and 0s - why do we care about this? Well because we can use this fact to do some checks
sum(metatable[,1] %in% colnames(intable)) # returns 46
length(metatable[,1] %in% colnames(intable)) # returns 48 - 48 entries, 46 of which are TRUE



## Now This big TRUE/FALSE list we can use as an index - where basically we use the TRUE and FALSE entries to indicate which of these we want to keep. 
metatablefilt = metatable[metatable[,1] %in% colnames(intable),]
intablefilt = intable[,colnames(intable) %in% metatablefilt[,1]]


## NEXT - lets try and filter our intable - its currently 20000 rows, but a lot of thos have a lot of 0s, so we want to subset this table to only include genes that have a minimum number of reads in a minimum number of samples
mingenecount = 4 # min number of reads in a gene for it to be significant
minsamples = round(ncol(intable)/2) # look from the inside out - take the number of columns in our intable, divide that by 2 (get half) then round that number (in case its a fraction) - so we want there to be the min number of reads in at least half of our samples



## Now that we have established our cutoffs - how can we use this to subset our data - so what are we doing again?
# keep ALL ROWS (genes) where NUMBER OF COLUMNS (samples) that have a COUNT HIGHER THAN MINGENECOUNT is GREATER THAN the MINSAMPLES
# so lets break this up
part1 = intablefilt >= mingenecount # run through every entry in our intable, and ask the question, is that entry greater than or equal to the mingenecount
part2 = rowSums(part1) # just like we added up out TRUEs before, we can do that again - this time, across each of the rows - to see how many TRUEs we have per row
part3 = part2 >= minsamples # how many of these genes have a count HIGHER than the minimum number of samples that we decided we needed to have
intablefilt2 = intablefilt[rowSums(intable >= mingenecount) >= minsamples,]



## FINALLY - let subset our data to only work with the T47D Cell Line - how would you do this?
# 1 - select METATABLE to only include ROWS that have "T47D" in the CELL COLUMN
metatablesubset = metatablefilt[metatablefilt[,"cell"] == "T47D",]
# 2 - select the INTABLE to only include COLUMNS that have the SAME LABEL as those in the FIRST COLUMN of our METATABLESUBSET
intablesubset = intablefilt2[,colnames(intablefilt2) %in% metatablesubset[,1]]


###########
## GRAPHICS
###########
## R has a lot of functionality with graphics. They have a lot of built in simple functions that will give you snapshots of what is going on. A widely held opinion though is the way you really get nice graphics - is by using GGPLOT2 - a WIDELY applicable graphic interface that allows you to make a lot of nice figures in a piecemeal fashion that allows for maximal customizability
# https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf
## Theres lots of thing you can do with ggplot2 - Im going to cover a couple here. But you can look up doing the google to see other things you can do
## Before we do anything - we have to download the ggplot2 package - this can be done either via CONDA or via the R COMMAND LINE
# install.packages("ggplot2")
# """ conda install r-ggplot2 """
library(ggplot2)


##### GGplot2 Scatterplot
## Starting with a scatterplot because its basic. 
## I picked a couple of potentially interesting genes (GOI = Genes Of Interest)
GOI = c("ESR1", "FOXA1", "GATA3", "TFAP2C", "PGR")

## Format our data - a good place to start for plotting is to have a table where samples are along the y-axis, and then each column is either a value to plot, or a piece of metadata to modify the plot (more on this later) - so lets make this table
## 1 - lets subset our data and pick out two genes - SELECT from INTABLESUBSET the rows of ESR1 and TFAP2C and ALL COLUMNS, and then TRANSPOSE the table
## LASTLY - AND THIS IS A WEIRD THING BUT NEEDS TO BE NOTED - you have to wrap the whole thing in the DATA.FRAME function to turn it into a data frame. This is a little outside the scope of a basic tutorial, but the two major tables are matrices and data.frames. most things are data.frames, matrices are numeric specific, and when we transposed our data - that automatically turned it into a matrix, so we have to turn it back into a data.frame
plottab1 = data.frame(t(intablesubset[c("ESR1", "TFAP2C"),]))
## This is a good start, but we are going to add a couple more columns of metadata on to the table - couple ways to do that
# 1 - cbind will push tables together - easy, but not the best
plottab1a = cbind.data.frame(plottab1, metatablesubset[,c("estrogen","dox")])
# 2 - the reason its not the best is because what if the order of the samples was different in our metatable and intable? We wouldnt know that, it would just mush it together, NOT GOOD, and an easy way to get errors
plottab1b = merge(x = plottab1, y = metatablesubset[,c("sample","estrogen","dox")], by.x = "row.names", by.y = "sample", sort = FALSE)
## Much safer! But notice that we not a column of samples - with the label as row.names - so we either want to get rid of that column, and add it back as rownames, or just rename and keep that column as labels, I think for readability and ease of access right now, lets just keep the column there, and rename it
colnames(plottab1b)[1] = "sample"

## Before anything - lets rename this back to plottab1
plottab1 = plottab1b




## ITERATION 1
#### So - ggplot2 works by establishing the AESTHETICs of the object, and then you build off of there. USE THE CHEATSHEET (everyone does...)
pout = ggplot(data = plottab1, mapping = aes(x = plottab1[,2], y = plottab1[,3]))
## From there - you built our the various pieces that you need for the plot - namely, we havent actually plotted anything - so we have to add POINTS to our plot to make it a scatter plot
pout = pout + geom_point()
print(pout)
## Now we can add labels to everything to make it clearer
# NOTE what I did for the y label. I forgot what was in the second column - but I know that string was already saved as the 2nd column name - so using this fact, I called the column names of the plot table, and then chose the second value. Not only is this nice because Im too lazy to look it up, but much more importantly - it makes this variable DYNAMIC - which means IT WILL CHANGE AUTOMATICALLY AS I ALTER THE UPSTREAM DATA - this is a fundamental tenant of good coding, and something to keep in mind. This parameter is SOFTCODED - and will adjust, as opposed to "ESR1" being hardcoded, and will not change even if i plot different data - which will obviously lead to mistakes!
pout = pout + labs(title = "scatter plot", x = "ESR1", y = colnames(plottab1)[3])
print(pout)

## ITERATION 2
## Now lets get a little fancier - I want to know which of these plots belong to which metadata category, add trendlines, alter shapes, label my points, etc.. so here is how you can do all of these things!
## Add labels that automatically shift away from the points
## Before we do anything - we have to download the ggplot2 package - this can be done either via CONDA or via the R COMMAND LINE
install.packages("ggrepel")
# """ conda install r-ggrepel """
library(ggrepel)
pout = ggplot(data = plottab1, mapping = aes(x = plottab1[,2], y = plottab1[,3], label = plottab1[,1]))
pout = pout + geom_point()
pout = pout + geom_text_repel(box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
print(pout)

## Add trendlines
## Define the function that will pull out the information in our trend line to plot on our figure
lm_eqn = function(df){
  colnames(df) = c("x", "y")
  m = lm(y ~ x, df);
  # eq = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
  #                  list(a = format(coef(m)[1], digits = 2),
  #                       b = format(coef(m)[2], digits = 2),
  #                       r2 = format(summary(m)$r.squared, digits = 3)))
  eq = substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), pvalue = format(summary(m)$coefficients[2,4], digits = 2)))
  as.character(as.expression(eq));
}

pout = ggplot(data = plottab1, mapping = aes(x = plottab1[,2], y = plottab1[,3], label = plottab1[,1]))
pout = pout + geom_point()
pout = pout + geom_smooth(data = plottab1, inherit.aes = FALSE, aes(x = plottab1[,2], y = plottab1[,3]), method = lm, color = "black")
pout = pout + geom_text(label = lm_eqn(plottab1[,c(2:3)]), parse = TRUE, x = 0.2*max(plottab1[,c(2)]), y = 1*max(plottab1[,c(3)]), size = 2)
print(pout)

## Add coloring to the points based off of the conditions
pout = ggplot(data = plottab1[,c(2,3,4)], aes(x = plottab1[,2], y = plottab1[,3], color = plottab1[,4], shape = plottab1[,5]))
pout = pout + geom_point()
pout = pout + labs(color = "Estrogen", shape = "DOX")
print(pout)


## THROW IT ALL TOGETHER:
## Add coloring to the points based off of the conditions
pout = ggplot(data = plottab1[,c(2,3,4)], aes(x = plottab1[,2], y = plottab1[,3], color = plottab1[,4], shape = plottab1[,5], label = plottab1[,1]))
pout = pout + geom_point()
pout = pout + geom_text_repel(box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"))
pout = pout + labs(title = "scatter plot", x = "ESR1", y = colnames(plottab1)[3], color = "Estrogen", shape = "DOX")
pout = pout + geom_smooth(data = plottab1, inherit.aes = FALSE, aes(x = plottab1[,2], y = plottab1[,3]), method = lm, color = "black")
pout = pout + geom_text(inherit.aes = FALSE, label = lm_eqn(plottab1[,c(2:3)]), parse = TRUE, x = 0.2*max(plottab1[,c(2)]), y = 1*max(plottab1[,c(3)]), size = 2)
print(pout)

## Saving a plot out - define the outfile, open it up, print the plot to it, and close the pdf
outfile = "/Users/mgc439/Code/codingclub/Rtutorial/outscatterplot.pdf"
pdf(outfile)
print(pout)
junk = dev.off()




##### BOXPLOTS
## So just like any ggplot graphics, you follow the same formula for a lot of the set up around the plotting - but then use a different call for the figure itself. So lets set up out data for boxplots
# First - we are going to make a new table, where we have our samples, our data, and then a column that encompasses all of our conditions concatenated together
# NEW TABLE = CBIND(mush together) the FIRST 3 COLUMNS of PLOTTAB and a 4th COLUMN, NAMED "condition" that is composed of a PASTED TOGETHER column of COLUMN 4 and 5 and SELECTING JUST THE ROWS THAT HAVE "noDOX" in the DOX condition (this will make it easier to work with later)
plottab2a = cbind.data.frame(plottab1[plottab1[,5]=="noDOX",c(1,2,3)], condition = paste(plottab1[plottab1[,5]=="noDOX",4], plottab1[plottab1[,5]=="noDOX",5], sep="_"))
# the melt function is a unique function - just do it and see what happens! Basically it "melts" the table so that the numerical values are all n a single column, and it creates a separate column to annotate the rows
# install.packages("reshape2")
# """ conda install r-reshape2 """
library("reshape2")
plottab2 = melt(plottab2a, id.vars = c("sample", "condition"))
plottab2[,3] = as.character(plottab2[,3])
plottab2[,2] = as.character(plottab2[,2])

## Now lets make our boxplot (less explanation this time - you can play with it on your own)
# as far as aes goes though - x is still along the xaxis (This time the category), y is along the y-axis - so the value, and z is the subcategory with the group
pout = ggplot(data = plottab2, aes(x = plottab2[,3], y = plottab2[,4], z = plottab2[,2], fill = plottab2[,2]))
pout = pout + geom_boxplot(na.rm=TRUE, outlier.shape = NA, position = position_dodge(width=0.9))
pout = pout + labs(x = "gene", y = "value", fill = "Estrogen", title = "Estrogen vs NoE2 for ESR1 and TFAP2C")
pout = pout + geom_point(data = plottab2, aes(x=plottab2[,3], y=plottab2[,4], fill = plottab2[,2]), 
                         position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 1)
print(pout)

## Now lets add signifigance levels (the things that everyone always wants)
install.packages("ggpubr")
# """ conda install r-ggpubr """
install.packages("ggsignif")
# """ conda install r-ggsignif """
library("ggpubr")
library("ggsignif")

## This is a useful function that will create all of the unique combinations of your list of items, that are 2 long - execute the command and see what happens
combtab = combn(unique(plottab2[,3]), 2, simplify=F)

pout = ggplot(data = plottab2, aes(x = plottab2[,3], y = plottab2[,4], z = plottab2[,2], fill = plottab2[,2]))
pout = pout + geom_boxplot(na.rm=TRUE, outlier.shape = NA, position = position_dodge(width=0.9))
pout = pout + labs(x = "gene", y = "value", fill = "Estrogen", title = "Estrogen vs NoE2 for ESR1 and TFAP2C")
pout = pout + geom_point(data = plottab2, aes(x=plottab2[,3], y=plottab2[,4], fill = plottab2[,2]), 
                          position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 1)
pout <- pout + stat_compare_means(aes(group = plottab2[,2]), method = "t.test", label = "p.format", label.y = max(plottab2[,4])*1.1, size = 3)
print(pout)


## Now lets add signifigance bars for comparisons between individual boxes
plottab3 = plottab2[plottab2[,3] == "TFAP2C",]
combtab = combn(unique(plottab3[,2]), 2, simplify=F)

pout = ggplot(data = plottab3, aes(x = plottab3[,2], y = plottab3[,4], fill = plottab3[,2]))
pout = pout + geom_boxplot(na.rm=TRUE, outlier.shape = NA, position = position_dodge(width=0.9))
pout = pout + labs(x = "gene", y = "value", fill = "Estrogen", title = "Estrogen vs NoE2 for ESR1 and TFAP2C")
pout = pout + geom_point(data = plottab3, aes(x=plottab3[,2], y=plottab3[,4], fill = plottab3[,2]), 
                         position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.12), size = 1)
pout <- pout + geom_signif(comparisons = combtab, step_increase = 0.1, y_position = c(max(plottab3[,4])*1.1), map_signif_level=FALSE, test = "t.test")
print(pout)


## Saving a plot out - define the outfile, open it up, print the plot to it, and close the pdf
outfile = "/Users/mgc439/Code/codingclub/Rtutorial/boxplot.pdf"
pdf(outfile)
print(pout)
junk = dev.off()


## Lastly - lets take on a new challenge - but one that everyone wants to know how to make - the HEATMAP!!!
dim(intablefilt2) # this is the count table we are going to use
dim(metatablefilt) # this is the metadata that we are going to put on our heatmap
## Now - heres the thing with super detailed ana beautiful heatmaps - they are HARD. more detail = more work, and although there are one liners out there to make heatmaps, if you really want to have full control over all the things you put into a heatmap - then you have to actually code each of those things you control. Along those lines, I think the best and most comprehensive heatmap package out there is ComplexHeatmap. People have different opinions - this is just mine, and I will show you how to use it
## Now lets add signifigance levels (the things that everyone always wants)

## NOTE - here that ComplexHeatmap - the package we are using, is from bioconductor - a large repository/organization of R packages specifically built for the biological sciences. Note that they use their own package manager - BiocManager, so we have to install that first
install.packages("BiocManager")
library(BiocManager)
install("ComplexHeatmap")
# """ conda install bioconductor-complexheatmap """
library("ComplexHeatmap")
# https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html - a VERY VERY COMPREHENSIVE AND GOOD GUIDE on how to use this.

## Now - we could just plot our table with a single call and be done with it - in fact, you can try that here:
# NOTE though that this will take a long time and not look very good.. why? run it and then read the reasons below

outhm = Heatmap(intablefilt2)
outhmfile = "/Users/mgc439/Code/codingclub/Rtutorial/outheatmap.pdf"
pdf(outhmfile)
print(outhm)
junk <- dev.off()

## Why did this take so long and look terrible
# 1 - you are putting in ALL of your data - most heatmaps will downsize their data to make it more manageable while capturing the same result
# 2 - the data wasnt scaled - ESPECIALLY FOR RNASEQ (but a lot of other applciations as well) you need to have your data scaled to put it in a heatmap - and by scaling, I mean going across each row of the data and rescoring it as a zscore - allowing for easy comparison across columns, but without having the raw numbers affect each other across rows.
# 3 - you didnt control the extraneous things like labels and such, which can be noisy or impossible to plot or can mess up an image in any number of ways

## So lets try this again - this time with several steps and parameters in the heatmap to make it look like something
# Step 1 - select the top 1000 genes based off of variance - common step and worth doing to downsize the data
# create new object - COUNTTAB_VARIANCE - that is a sorted list of the VARIANCE per row, VARIANCE per ROW is determined using the inner code segment
# apply(intablefilt2, 1, var) - APPLY the function VAR (variance) to each ROW (1) of the object INTABLEFILT2
counttab_variance = sort(apply(intablefilt2,1,var), decreasing=TRUE)
# starting with INTABLEFILT2, select only the NAMES (genes) of the top 1:1000 VARIANCES and all columns
hmplottab = intablefilt2[names(counttab_variance[1:1000]),]

# Scale each row via zscore
# Define a quick function that when applied to a vector of numbers - will create the zscore
zscore <- function(x) {(x-mean(x))/sd(x)}
# Then apply that function (like before) to each row of our data - NOTE that this function when applied to a whole table, but returns a VECTOR of data, builds out on a column by column basis, essentially transposing our matrix, so we have to transpose it back
# NOTE - that we also have to assign it as a matrix for later use, sorry!
hmplottabscaled = as.matrix(t(apply(hmplottab, 1, zscore)))
#as.matrix(t(apply(counttab, 1, function(x) zscore(x))))

## Now that we have our data the way we like it, lets plot it - but with some parameters
outhm = Heatmap(hmplottabscaled,
                row_title = "Genes",                                       ## Name the rows
                column_title = "Samples",                                  ## Name the columns
                
                cluster_columns = TRUE,                         ## Cluster the columns or leave as is
                cluster_rows = TRUE,                            ## Cluster the rows or leave as is
                
                show_column_names = TRUE,                                  ## Show the Column Names
                column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
                show_row_names = FALSE,                                    ## Show the row names
                row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
                row_names_gp = gpar(fontsize=6),
                
                show_row_dend = TRUE,                                       ## Show the dendrogram on the rows
                show_column_dend = TRUE,                                    ## Show the dendrogram on the columns
                
                heatmap_legend_param = list(title = "Zscore",
                                            legend_height = unit(2, "cm"),
                                            title_gp = gpar(fontsize = 8, fontface = "bold"))
                )

outhmfile = "/Users/mgc439/Code/codingclub/Rtutorial/outheatmap.pdf"
pdf(outhmfile)
print(outhm)
junk <- dev.off()



## That's more like a recognizable heatmap. now lets work on getting those annotations along the top of the heatmap to incorporate our metadata. This is where we start to get pretty tricky, and bugs will happen. We first need our metadata to match up with our mapping data - namely, we need to have the sample names that were the column names of our mapping data - to be row names for our metadata, and we need to have all the same names in each of those lists
# remind ourselves what our metadata looks like:
head(metatablefilt)
# make our first column of names be rownames
annotationtab = metatablefilt[,2:ncol(metatablefilt)]
rownames(annotationtab) = metatablefilt[,1]

# Next - our heatmap function is going to cluster our data - so we need to make sure that the metadata and mapping data are in the same place, so that when our function clusters the mapping data, it can take that cluster order - and then apply it to our metadata. IF YOUR ANNOTATION AND MAPPING DATA ARE NOT THE SAME ORDER TO START - THEN YOU WILL GET MISORDER ANNOTATIONS
# so lets do a simple check
all.equal(rownames(annotationtab), colnames(hmplottabscaled))

## Now thats taken care of - lets build our annotation
hatop = HeatmapAnnotation(df = annotationtab)

outhm = Heatmap(hmplottabscaled,
                row_title = "Genes",                                       ## Name the rows
                column_title = "Samples",                                  ## Name the columns
                
                cluster_columns = TRUE,                         ## Cluster the columns or leave as is
                cluster_rows = TRUE,                            ## Cluster the rows or leave as is
                
                show_column_names = TRUE,                                  ## Show the Column Names
                column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
                show_row_names = FALSE,                                    ## Show the row names
                row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
                row_names_gp = gpar(fontsize=6),
                
                show_row_dend = TRUE,                                       ## Show the dendrogram on the rows
                show_column_dend = TRUE,                                    ## Show the dendrogram on the columns
                
                heatmap_legend_param = list(title = "Zscore",
                                            legend_height = unit(2, "cm"),
                                            title_gp = gpar(fontsize = 8, fontface = "bold")),
                
                top_annotation = hatop
)

outhmfile = "/Users/mgc439/Code/codingclub/Rtutorial/outheatmap.pdf"
pdf(outhmfile)
print(outhm)
junk <- dev.off()













## WRITE HEATMAP FUNCTION - fully finished and functional - This is a highly customized and functionalized heatmap function that allows you to do a number of things by just changing the input parameters. This includes column AND row annotations, turning on and off the clstering, subsetting the number of top genes to use, etc. NOTE that it has a bunch of hidden automatic things to that will adjust the heatmap to fit funky names and data. It's not perfect (functions never are and they will break with something) - but it allows you to do a lot of heatmaps quickly by just changing the parameeters
## Load in Libraries
install.packages("viridis")
packagelist = c("ggplot2", "reshape2", "grid", "gridExtra", "scales", "ggrepel", "circlize", "viridis", "RColorBrewer", "tools", "ggpubr", "ggsignif")
junk <- lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

#counttab ## the counttable
#subsetnum = FALSE ## set to a number if you want to select by the N most varied genes (good for subsetting)
#metatable ## The metatable
#annotationlist ## The annotation object created beforehand (NEED TO AUTOMATE THIS MORE)
#colclusterparam = FALSE ## if not FALSE - will cluster the columns
#rowclusterparam = FALSE ## if not FALSE - will cluster the rows
#pdfoutfile ## the path to the outfile that the pdf will be saved to
create_heatmap <- function(counttab, samplesample = FALSE, subsetnum = FALSE, colmetatable = NULL, colannotationlist = NULL, 
                           rowmetatable = NULL, rowannotationlist = NULL, colclusterparam = FALSE, rowclusterparam = FALSE, pdfoutfile) {
  
  ## Filter to the N most varied genes if applicable
  if (subsetnum !=FALSE) {
    counttab_variance = sort(apply(counttab,1,var), decreasing=TRUE)
    counttab = counttab[names(counttab_variance[1:subsetnum]),]
  }

  ## FAILSAFE - if there are any columns that sneak in somehow which have a variance of 0 - they have to be excluded for plotting
  nearZeroVarcols <- which(sapply(counttab, var) == 0)
  if (length(nearZeroVarcols) > 0) {
    counttab = counttab[,!colnames(counttab) %in% names(nearZeroVarcols)]
    print(paste("excluding column ", names(nearZeroVarcols), " due to zero variance", sep=""))}

  
  ## Calc. spearman correlation and use values for column clustering before any other alterations
  #if (colclusterparam != FALSE) {
    ## FAILSAFE - if there are any columns that sneak in somehow which have a variance of 0 - they have to be excluded for plotting
    nearZeroVarcols <- which(sapply(counttab, var) == 0)
    if (length(nearZeroVarcols) > 0) {
      counttab = counttab[,!colnames(counttab) %in% names(nearZeroVarcols)]
      print(paste("excluding column ", names(nearZeroVarcols), " due to zero variance", sep=""))}
    cordata <- cor(counttab, method="spearman")
    coldistance = dist(t(as.matrix(na.omit(cordata))), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
    #colclusterparam = colcluster
  #}

  ## Zscore out counttable, or turn into spearman correlation if doing sample-sample comparison
  maptab = as.matrix(t(apply(counttab, 1, function(x) zscore(x))))

  #if (rowclusterparam != FALSE) {
    rowdistance = dist(maptab, method = "euclidean")
    rowcluster = hclust(rowdistance, method = "ward.D2")
    #rowclusterparam = rowcluster
  #}

  if (samplesample !=FALSE) {
    maptab <- cor(counttab, method="spearman")
    rowdistance = dist(as.matrix(cordata), method = "euclidean")
    rowcluster = hclust(rowdistance, method = "ward.D2")
    coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
  }
  
  if (rowclusterparam != FALSE) {rowclusterparam = rowcluster}
  if (colclusterparam != FALSE) {colclusterparam = colcluster}

  ## Build Annotations from our metatable
  hatop = NULL
  if (!is.null(colannotationlist) & !is.null(colmetatable)) {
    temp1 <- vector("list", length(colannotationlist))
    names(temp1) = names(colannotationlist)
    annotlegendlist = lapply(temp1, function(x) x[[1]] = list(title_gp=gpar(fontsize=5, fontface="bold"), labels_gp=gpar(fontsize=4)))
    ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
    showlegendparam = unname(unlist(lapply(colannotationlist, function(x) {
      numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
      is.null(numterms) || numterms <= 10})))
    hatop = HeatmapAnnotation(df = colmetatable,
                              col = colannotationlist,
                              na_col = "white",
                              show_annotation_name = TRUE,
                              annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                              annotation_name_side = "left",
                              height = unit(5, "cm"),
                              show_legend = showlegendparam,
                              annotation_legend_param = annotlegendlist)
  }

  ## Defune the side annotation if data is supplied
  if (!is.null(rowannotationlist) & !is.null(rowmetatable)) {
    
    ## Define parameters for each of the labels on the annotation bars
    temp1 <- vector("list", length(rowannotationlist))
    names(temp1) = names(rowannotationlist)
    annotlegendlist = lapply(temp1, function(x) x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))
    
    ## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
    showlegendparam = unname(unlist(lapply(rowannotationlist, function(x) {
      numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
      is.null(numterms) || numterms <= 10})))
    
    ## Look for any empty annotations - fill them with white, and later, make sure to hide their legend
    emptyannots = names(sapply(rowannotationlist, length)[sapply(rowannotationlist, length)==0])
    if (length(emptyannots) > 0){
      for (i in 1:length(emptyannots)) {
        temp1 = "white"
        names(temp1) = emptyannots[i]
        rowannotationlist[[emptyannots[i]]] = temp1
      }
      showlegendparam[which(names(rowannotationlist) %in% emptyannots)] = FALSE
    }
    haside = rowAnnotation(df = rowmetatable,
                           col = rowannotationlist,
                           na_col = "white",
                           gp = gpar(fontsize = 0.5),
                           show_annotation_name=TRUE,
                           annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                           annotation_name_side = "top",
                           show_legend = showlegendparam,
                           annotation_legend_param = annotlegendlist)
  }

  ## Define the Heatmap
  if (samplesample == FALSE) {heatmapcolorparam = colorRamp2(c(-3,0,3), c("blue", "white", "red"))} else{heatmapcolorparam = colorRamp2(c((min(maptab)),1), c("white", "red"))}
  ht1 = Heatmap(maptab,
                col = heatmapcolorparam,    ## Define the color scale for the heatmap
                row_title = "Genes",                                       ## Name the rows
                column_title = "Samples",                                  ## Name the columns
                
                cluster_columns = colclusterparam,                         ## Cluster the columns or leave as is
                cluster_rows = rowclusterparam,                            ## Cluster the rows or leave as is
                
                show_column_names = TRUE,                                  ## Show the Column Names
                column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
                show_row_names = nrow(maptab) <=200,                                    ## Show the row names
                row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
                row_names_gp = gpar(fontsize=6),
                
                show_row_dend = nrow(maptab) <=2000,                                     ## Show the dendrogram on the rows
                show_column_dend = TRUE,                                   ## Show the dendrogram on the columns
                
                heatmap_legend_param = list(title = ifelse(samplesample==FALSE, "Zscore", "Spearman\nCorrelation"),
                                            legend_height = unit(2, "cm"),
                                            title_gp = gpar(fontsize = 8, fontface = "bold")),
                top_annotation = hatop
                
  )

  ## Plot out the heatmap
  pdf(file = pdfoutfile, width=11,height=8.5)
  if (!is.null(rowannotationlist) & !is.null(rowmetatable)) {draw(ht1 + haside, annotation_legend_side = "bottom", padding = unit(c(5, 20, 5, 5), "mm"))} else {draw(ht1, annotation_legend_side = "bottom", padding = unit(c(5, 20, 5, 5), "mm"))}
  junk <- dev.off()
  
}



## Input the metatable and this will build an annotation, if you want to preset the colors - then input them in proper list form, either as named list of named colors, or with color brewer
# this is an example of a customcolorlist that you can input - where you either name the colors you want for each label in a metadata column, or create a scale of colors using colorRamp2 colorRamp2(EACH VALUE OF THE SCALE YOU WANT TO CREATE (ex: c(0,1,2)) , brewer.pal(NUMBER OF BREAKS, "COLOR OBJECT")
#customcolorlist = list(platelet_group = c("Hyper" = "red", "Hypo" = "blue"), 
#                       PF4 = colorRamp2(seq(min(metatable[,3], na.rm = TRUE), max(metatable[,3], na.rm = TRUE), length = 3), brewer.pal(3, "Purples"))
annotationlist_builder <- function(metatable, customcolorlist = NULL) {
  annotlist <- list()
  colorlist = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  for (colnum in 1:ncol(metatable)) {
    # if (!is.null(customdiscretecolorlist)) {
    #   collist = customdiscretecolorlist[colnames(metatable[,colnum])]}
    #   break
    # }
    if (!is.null(customcolorlist) && colnames(metatable[,colnum,drop=FALSE]) %in% names(customcolorlist)) {
      annotlist[colnames(metatable[,colnum,drop=FALSE])] = customcolorlist[colnames(metatable[,colnum,drop=FALSE])]}
    
    if (!is.null(customcolorlist) && !colnames(metatable[,colnum,drop=FALSE]) %in% names(customcolorlist) | is.null(customcolorlist)) {
      if (is.character(metatable[,colnum])) {
        #annotlist[[colnum]] = sample(colorlist, length(unique(metatable[,colnum])))
        annotlist[[colnum]] = sample(colorlist, length(na.omit(unique(metatable[,colnum]))))
        #names(annotlist[[colnum]]) = unique(metatable[,colnum])
        names(annotlist[[colnum]]) = na.omit(unique(metatable[,colnum]))
        names(annotlist)[[colnum]] = colnames(metatable)[colnum]
      }
      else {collist = colorRamp2(seq(min(metatable[,colnum], na.rm = TRUE), max(metatable[,colnum], na.rm = TRUE), length = 3), brewer.pal(3, "Purples"))
      annotlist[[colnum]] = collist
      names(annotlist)[colnum] = colnames(metatable)[colnum]
      }
    }
  } 
  return(annotlist)
}

## Now run it!
annotationlist1 = annotationlist_builder(annotationtab)
outhmfiletest = "/Users/mgc439/Code/codingclub/Rtutorial/outheatmap_test.pdf"
create_heatmap(intablefilt2, samplesample = FALSE, subsetnum = 1000, colmetatable = annotationtab, colannotationlist = annotationlist1,
               colclusterparam = TRUE, rowclusterparam = TRUE, pdfoutfile = outhmfiletest)

outhmfileSStest = "/Users/mgc439/Code/codingclub/Rtutorial/outheatmap_SS_test.pdf"
create_heatmap(counttab = intablefilt2, samplesample = TRUE, subsetnum = 1000, colmetatable = annotationtab, colannotationlist = annotationlist1,
               colclusterparam = TRUE, rowclusterparam = TRUE, pdfoutfile = outhmfileSStest)


