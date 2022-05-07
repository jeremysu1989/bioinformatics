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



