## R basic
### Everything in R is an object
#### R has 6 basic data types. (In addition to the five listed below, there is also raw which will not be discussed in this workshop.)

1 character

2 numeric (real or decimal)

3 integer

4 logical

5 complex

#### R provides many functions to examine features of vectors and other objects, for example

class() - what kind of object is it (high-level)?

typeof() - what is the object’s data type (low-level)?

length() - how long is it? What about two dimensional objects?

attributes() - does it have any metadata?

#### R has many data structures. These include

atomic vector

list

matrix

data frame

factors

#### Objects can have attributes. Attributes are part of the object. These include:

names

dimnames

dim

class

attributes (contain metadata)

#### In R matrices are an extension of the numeric or character vectors. They are not a separate type of object but simply an atomic vector with dimensions; the number of rows and columns. As with atomic vectors, the elements of a matrix must be of the same data type.

#### In R lists act as containers. Unlike atomic vectors, the contents of a list are not restricted to a single mode and can encompass any mixture of data types. Lists are sometimes called generic vectors, because the elements of a list can by of any type of R object, even lists containing further lists. This property makes them fundamentally different from atomic vectors.

#### A data frame is a very important data type in R. It’s pretty much the de facto data structure for most tabular data and what we use for statistics.  A data frame is a special type of list where every element of the list has same length (i.e. data frame is a “rectangular” list).

![image](https://user-images.githubusercontent.com/104820908/167233994-bfe7bbee-0ae4-4828-8baf-052e62c26a54.png)

#### Key Points
R’s basic data types are character, numeric, integer, complex, and logical.

R’s basic data structures include the vector, list, matrix, data frame, and factors. Some of these structures require that all members be of the same data type (e.g. vectors, matrices) while others permit multiple data types (e.g. lists, data frames).

Objects may have attributes, such as name, dimension, and class.


## reshape data

tidyr

## dataframe slicing and dicing
rm(list = ls())

x <- sample(1:50, 300, replace = TRUE)

y <- 3.2*x + rnorm(300,0,40)

d <- data.frame(x=x, y=y)

summary(d$x)

summary(d$y)

d$x >=48

d[d$x >= 48,]

d[d$x >=45 & d$y >=200,]

d[d$x >=45 & d$y >=200,c("x")]

d[d$x >=45 & d$y >=200,c("y")]

d$x[d$y > 200]

summary(d$x[d$y > 200])

which(d$y > 200)

which(d$x >=45 & d$y >=200)

which.min(d$x)

d[which.min(d$x),]

which.max(d$y)

d[which.max(d$y),]

subset(d,d$x <= 2)

## Binnind data with cut() and Bar plots with ggplot2
library("ggplot2")

d$x.binned <- cut(d$x,5)

d$x.binned

table(d$x)

table(d$x.binned)

ggplot(d) + geom_bar(aes(x=x))

ggplot(d) + geom_bar(aes(x=x.binned, fill=x.binned))

## Matching vectors and merging dataframes

c(3, 4, -1) %in% c(1, 3, 4, 8)

match(c("A", "C", "E", "A"), c("A", "B", "A", "E"))

merge(x, y, by = intersect(names(x), names(y)),
      by.x = by, by.y = by, all = FALSE, all.x = all, all.y = all,
      sort = TRUE, suffixes = c(".x",".y"), no.dups = TRUE,
      incomparables = NULL, ...)

## ggplot2 Facet

https://ggplot2-book.org/facet.html

http://www.sthda.com/english/wiki/ggplot2-facet-split-a-plot-into-a-matrix-of-panels

## List
adh <- list(chr="2L", start=14615555L, end=14618902L, name="Adh")

adh$chr

adh[1]

adh[[1]]

adh['chr']

adh[['chr']]

adh$id <- "FBgn0000055"

adh$id2 <- c("FBgn0000055","FBgn0000077")

adh$id2 <- NULL

unlist(adh)

## lapply() and sapply()
#### calculate the mean value of each vector in list with loop
ll <- list(a=rnorm(6, mean=1), b=rnorm(6, mean=4), c=rnorm(6, mean=6))

ll_means <- numeric(length(ll))

for (i in seq_along(ll)){
  ll_means[i] <- mean(ll[[i]])
}

typeof(ll_means)       #"double"
#### calculate the mean value of each vector in list with lapply()
ll_means <- lapply(ll,mean)

typeof(ll_means)       #"list"

#### specify additional argument
ll_means <- lapply(ll, mean, na.rm=TRUE)

#### Or use anonymous funciton
ll_means <- laply(ll, function(x) mean(x, na.rm=TRUE))

#### write a more verbose function
      meanRemoveNAVerbose <- function(x, warn=TRUE) {
            #A function that removes missing values when calculating the mean     
            #and warns us about it.
            if (any(is.na(x)) && warn) {
                  warning("removing some missing values!")
            }
            mean(x, na.rm=TRUE)
      }

#### Working with the Split-Apply-Combine Pattern

There are a few other useful tricks to know about the split-apply-combine pattern
built from split(), lapply(), and do.call() with rbind()
#### Split
d_split <- split(d$x,d$x.binned)

str(d_split)
#### Apply
lapply(d_split,summary)
#### Combine
df1 <- do.call(cbind,lapply(d_split,summary))

df2 <- do.call(rbind,lapply(d_split,summary))

## Exploring Dataframe with pdlyr
dplyr has five basic functions for manipulating dataframes: arrange(), filter(),
mutate(), select(), and summarize()

## working with string
grep, grepl, regexpr, gregexpr, regexec and gregexec search for matches to argument pattern within each element of a character vector: they differ in the format of and amount of detail in the results.

#### sub and gsub perform replacement of the first and all matches respectively.

sub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE,
    fixed = FALSE, useBytes = FALSE)
    
      str <- "Now is the time      "
      sub(" +$", "", str)     
      
      sub("gene=(\\w+)", "\\1", "gene=LEAFY", perl=TRUE)
      
      ## capitalizing
      txt <- "a test of capitalizing"
      gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", txt, perl=TRUE)
      gsub("\\b(\\w)",    "\\U\\1",       txt, perl=TRUE)
#### paste    
paste (..., sep = " ", collapse = NULL, recycle0 = FALSE)

paste0(...,            collapse = NULL, recycle0 = FALSE)

    paste0(1:12, c("st", "nd", "rd", rep("th", 9)))
    paste(1:12, c("st", "nd", "rd", rep("th", 9)))
    
 #### strspilt  
strsplit(x, split, fixed = FALSE, perl = FALSE, useBytes = FALSE)

      leafy <- "gene=LEAFY;locus=2159208;gene_model=AT5G61850.1"
      strsplit(leafy, ";")

## Developing Workflows with R Scripts
#### Control Flow: if, for, and while

The basic syntax of if, for, and while are:

      if (x == some_value) {
            # do some stuff in here
      } else {
            # else is optional      
      }
      for (element in some_vector) {
       # iteration happens here
      }
      while (something_is_true) {
            # do some stuff
      }
 #### Working with R Scripts     
You can run R scripts from R using the function source(). Also, it’s a good idea to
indicate (either in comments or a README file) which directory the user should set
as their working directory

       source("my_analysis.R")

#### Workflows for Loading and Combining Multiple Files
In bioinformatics projects, sometimes loading data into R can be half the battle.
Loading data is nontrivial whenever data resides in large files and when it’s scattered
across many files. These workflows tie together many of the R
tools we’ve learned so far: lapply(), do.call(), rbind(), and string functions like sub().

The first step to loading this data into R is to programmatically access these files from
R. If your project uses a well-organized directory structure and consistent filenames, 
this is easy. R’s function list.files() lists all files in a specific directory. 

With list.files() programmatically listing our files, it’s then easy to lapply() a
function that loads each file in. 
      loadFile <- function(x) {
            # read in a BED file, extract the chromosome name from the file,
            # and add it as a column
            df <- read.delim(x, header=FALSE, col.names=bedcols)
            df$chr_name <- sub("hotspots_([^\\.]+)\\.bed", "\\1", basename(x))
            df$file <- x
            df
      }
            
Just as there’s more than one way to pet a cat, there are many ways to bulk load data into R.                 
Processing files this way is very convenient for large files, and is even more powerful
because it’s possible to parallelize data processing by simply replacing lapply() with
mclapply().

#### Exporting Data

      write.table(x, file = "", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
 
save writes an external representation of R objects to the specified file. The objects can be read back from the file at a later date by using the function load or attach (or data in some cases).

save.image() is just a short-cut for ‘save my current workspace’, i.e., save(list = ls(all.names = TRUE), file = ".RData", envir = .GlobalEnv). It is also what happens with q("yes")

      save(..., list = character(),
            file = stop("'file' must be specified"),
            ascii = FALSE, version = NULL, envir = parent.frame(),
            compress = isTRUE(!ascii), compression_level,
            eval.promises = TRUE, precheck = TRUE)

      save.image(file = ".RData", version = NULL, ascii = FALSE,
           compress = !ascii, safe = TRUE)
