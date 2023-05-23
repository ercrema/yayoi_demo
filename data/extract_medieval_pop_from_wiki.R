library(rvest)

# Set the URL for the wiki page
url <- "https://ja.wikipedia.org/wiki/%E8%BF%91%E4%BB%A3%E4%BB%A5%E5%89%8D%E3%81%AE%E6%97%A5%E6%9C%AC%E3%81%AE%E4%BA%BA%E5%8F%A3%E7%B5%B1%E8%A8%88"

# Read the HTML code from the website
webpage <- read_html(url)

# Extract the tables from the webpage
tables <- webpage |> html_table(fill = TRUE)
estimates=tables[[5]]
sinK <- apply(estimates[,2:ncol(estimates)],2,remK)
sinK=apply(sinK,2,unlist)
estimates[,2:ncol(estimates)]=apply(sinK,2,as.numeric)
