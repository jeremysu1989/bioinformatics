# this script is used for scrape hyper-links from the website as follow
# https://www.geeksforgeeks.org/cross-validation-in-r-programming/

rm(list=ls())
if(!require(pacman)){
  install.packages('pacman')
  require(pacman)
}
pacman::p_load(xml2,rvest,dplyr,stringr,magrittr,curl, tidyverse)

# Start by reading a main test with read_html():
webpage <- read_html("https://www.geeksforgeeks.org/cross-validation-in-r-programming/")

# try to get the html-links related titles
main_text <- webpage %>% html_nodes("li a") %>% html_text2()
main_text

# try to get the hyper link in the whole page
wholelink <- webpage %>% html_nodes("li a") %>% html_attrs()
wholelink

# try to get the text of the table of contents
tabletext <- webpage %>% html_nodes("#home-page > div > div.sideBar > div > ul") %>% html_text2()
tabletext
tabletext %>% str_split("\n")

# try to get the hyper link of the table of contents
tablelink <- webpage %>% html_nodes("#home-page > div > div.sideBar > div > ul") %>% html_nodes("a") %>% html_attrs()
tablelink %>% unlist %>% length()

