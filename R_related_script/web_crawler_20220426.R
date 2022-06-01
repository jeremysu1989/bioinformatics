# get the poster information from Grail website
# install and load packages
rm(list=ls())
if(!require(pacman)){
  install.packages('pacman')
  require(pacman)
}
pacman::p_load(xml2,rvest,dplyr,stringr,magrittr,curl)


##### - Get the HTML for the web page that you want to scrape
##### - Decide what part of the page you want to read and find out what HTML/CSS you need to select it
##### - Select the HTML and analyze it in the way you need

# collect the name of the fiction and download times
# generate a new empty dataframe
df=as.data.frame(matrix(nrow = 0,ncol = 0))


## prepare the session
url <- "https://movie.douban.com/"
session <- session(url)

## prepare the form
forms <- html_form(read_html(url))
forms
form <- forms[[1]]
form

# add the search text
filled_form <- html_form_set(form, search_text = "ÃÀÈËÓã")
filled_form
## submit the session to web server
session2 <- html_form_submit(filled_form)

## test the setting
iconv(URLdecode(session2$url),"UTF8")

## get the image
desired_page <- session2 %>% read_html()
img_link <- desired_page %>% html_nodes("div.item-root  a")
img_link
