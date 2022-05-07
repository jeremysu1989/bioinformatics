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














