# get the poster information from AACR
# install and load packages
rm(list=ls())
if(!require(pacman)){
  install.packages('pacman')
  require(pacman)
}
pacman::p_load(xml2,rvest,dplyr,stringr,magrittr)

##### - Get the HTML for the web page that you want to scrape
##### - Decide what part of the page you want to read and find out what HTML/CSS you need to select it
##### - Select the HTML and analyze it in the way you need

# collect the name of the fiction and download times
# generate a new empty dataframe
df=as.data.frame(matrix(nrow = 0,ncol = 0))

# use for loop to access each webpege
for(i in (2:10)){url <- paste("https://www.txt80.com/hot/index_",i,".html",sep = "")
  web <- read_html(url)
# use each command to fetch the target information
# use *magrittr* to link different commands 
  title <- web %>% html_nodes("h4 a") %>% html_text() %>% str_replace_all("全本TXT小说下载","")

  downld <- web %>% html_nodes("span b") %>% html_text()
  downld <- downld[2:length(downld)] %>% str_replace_all("人下载","")

  cls <- web %>% html_nodes("b a") %>% html_text()
  fic_type <- cls[seq(1,length(cls),by =2)]
  fic_author <- cls[seq(2,length(cls),by =2)]

  pub_time <- web %>% html_nodes("p font") %>% html_text()
  pub_time <- pub_time[seq(2,length(pub_time),by =2)] %>% str_replace_all("发布时间：","")

  intro <- web %>% html_nodes("p") %>% html_text()
  intro1 <- intro[seq(1,length(intro),by =3)]
  # to split the vector, get list, use unlist to generate vector
  fic_size <- intro1[3:length(intro1)-1] %>% strsplit(split = "：") %>% unlist
  fic_size <- fic_size[seq(5,length(fic_size),by=5)]
  intro3 <- intro[seq(3,length(intro),by =3)]
  fic_intro <-intro3[1:length(intro3)-1]

#### combine these information
  fic_df <- as.data.frame(cbind(title,downld,fic_type,fic_author,pub_time,fic_size,fic_intro))
  fic_df$webpage = i
#### combind different fic_df in each iteration  
  df = combine(df,fic_df)
}
dim(df)
