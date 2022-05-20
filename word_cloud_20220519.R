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

###############################################################################
# the following code is prepared based on https://www.jianshu.com/p/9130e9349e20
wordcloud2(data = d)
wordcloud2(data = d, color = "random-light", backgroundColor = "white",shape = 'pentagon')
# change the shape of the word cloud
# wordcloud2(data = d, figPath = 'cancer.jpg', size = 1.5)
# change the shape to letter
# letterCloud(data = d, word = "R", size = 2)

###############################################################################
# for Chinese word, please confirm the fileEncoding = ‘utf8’, in case the appearance of unreadable or messy code
wordcloud2(data, size = 1, shape='cardioid',color = 'random-dark', backgroundColor = "pink",fontFamily = "微软雅黑")

##############################################################################
# for more detail or verbose examples, please visit:
# https://zhuanlan.zhihu.com/p/499935064
# https://bookdown.org/Maxine/tidy-text-mining/
# https://scholar.harvard.edu/files/jbenchimol/files/text-mining-methodologies.pdf
# https://blog.csdn.net/gdyflxw/article/details/53535227
# https://www.jianshu.com/p/806cdccd7a9a
