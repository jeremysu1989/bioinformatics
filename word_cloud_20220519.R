# word cloud
# this script is used for word frequency analysis in English.
# the first part is based on http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know/

rm(list=ls())
if(!require(pacman)){
  install.packages('pacman')
  require(pacman)
}

# install the related packages
p_load(tm,SnowballC,wordcloud,RColorBrewer)


#########Text mining
# set the working directory
getwd()
setwd("C:/Users/jeremy.su/Documents")
# load the text
text <- readLines("abstract.txt", encoding = 'ANSI')
# load the data as a corpus
docs <- Corpus(VectorSource(text))
# inspect the content of the document
inspect(docs)
# test transformation
# transformation is performed using tm_map() function to replace special characters from the text
# replacing "/" "@" and "|" with space
toSpace <- content_transformer(function(x, pattern) gsub(pattern, "",x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "|")
# cleaning the text
# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove your own stop word
# specify your stopwords as a character vector
docs <- tm_map(docs, removeWords, c("samples", "using", "act", "may", "compared", "however"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)

# build a term-document matrix
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 20)

# generate the word cloud
set.seed(1231)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=100, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))

