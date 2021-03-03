params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

knitr::opts_chunk$set(echo = TRUE)

#knitr::opts_knit$set(root.dir = "RU_tidyverse/inst/extdata/")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Making R easier with Tidyverse

---
"    
  )
  
}



## -----------------------------------------------------------------------------

load(file='data/my_tidy.Rdata')


## -----------------------------------------------------------------------------
head(df1)


## -----------------------------------------------------------------------------
head(df2)


## -----------------------------------------------------------------------------
head(df3a)


## -----------------------------------------------------------------------------
head(df3b)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Reading and Tibbles

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Reading and Tibbles

---
"    
  )
  
}



## -----------------------------------------------------------------------------
library(tidyverse)


## -----------------------------------------------------------------------------
untidy_counts_base <- read.csv("data/hemato_rnaseq_counts.csv")

untidy_counts_base


## -----------------------------------------------------------------------------

read_csv("data/hemato_rnaseq_counts.csv")



## -----------------------------------------------------------------------------
untidy_counts <- read_csv("data/hemato_rnaseq_counts.csv", col_types = cols(
    ENTREZ = col_character(),
    CD34_1 = col_integer(),
    ORTHO_1 = col_integer(),
    CD34_2 = col_integer(),
    ORTHO_2 = col_integer()
  ))
untidy_counts


## -----------------------------------------------------------------------------
untidy_counts[1,]

## -----------------------------------------------------------------------------
untidy_counts[,1]


## -----------------------------------------------------------------------------

untidy_counts[1]


## -----------------------------------------------------------------------------
untidy_counts[[1]]


## -----------------------------------------------------------------------------
untidy_counts$ENTREZ


## -----------------------------------------------------------------------------

as_tibble(untidy_counts_base)


## -----------------------------------------------------------------------------

untidy_counts_base <- as_tibble(untidy_counts_base)
untidy_counts_base <- mutate_at(untidy_counts_base, vars(ENTREZ), as.character)
untidy_counts_base


## -----------------------------------------------------------------------------

as.data.frame(untidy_counts_base) %>% head(n=12)


## -----------------------------------------------------------------------------
# Lets load in some packages
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

hg19_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

keys <- hg19_genes$gene_id

symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = keys, columns = c("SYMBOL"), keytype = "ENTREZID")



## -----------------------------------------------------------------------------

counts_metadata <- tibble(ID = symbols$ENTREZID, SYMBOL = symbols$SYMBOL, LENGTH = lengths(hg19_genes))

counts_metadata


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Tidying up your data

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Tidying up your data

---
"    
  )
  
}



## -----------------------------------------------------------------------------
untidy_counts


## -----------------------------------------------------------------------------
untidy_counts


## -----------------------------------------------------------------------------

tidier_counts <- pivot_longer(untidy_counts, names_to = "Sample", values_to = "counts", cols = c(-ENTREZ))
tidier_counts



## -----------------------------------------------------------------------------
pivot_wider(tidier_counts, names_from = c(Sample), values_from = counts)


## -----------------------------------------------------------------------------
tidier_counts


## -----------------------------------------------------------------------------
tidier_counts


## -----------------------------------------------------------------------------
tidy_counts <- separate(tidier_counts, Sample, sep = "_", into=c("CellType", "Rep"), remove=TRUE)
tidy_counts


## -----------------------------------------------------------------------------
unite(tidy_counts, Sample, CellType, Rep, remove=FALSE)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Piping with Magrittr

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Piping with Magrittr

---
"    
  )
  
}



## -----------------------------------------------------------------------------
tidy_counts <- separate(pivot_longer(untidy_counts, names_to = "Sample", values_to = "counts", cols = c(-ENTREZ)), Sample, sep = "_", into = c("CellType","Rep"), remove=FALSE)


## -----------------------------------------------------------------------------
tidier_counts <- pivot_longer(untidy_counts, names_to = "Sample", values_to = "counts", cols = c(-ENTREZ))
tidy_counts <- separate(tidier_counts, Sample, sep = "_", into = c("CellType","Rep"), remove=FALSE)



## -----------------------------------------------------------------------------

tidy_counts <- untidy_counts %>% 
  gather(key=Sample, value=counts, -ENTREZ) %>% 
  separate(Sample, sep = "_", into = c("CellType","Rep"), remove=FALSE)
tidy_counts


## -----------------------------------------------------------------------------

tidier_counts %>% .[.$counts > 0,]



## -----------------------------------------------------------------------------
library(magrittr)
tidier_counts %<>% .[.$counts > 0,]
tidier_counts


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Joining tibbles together

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Joining tibbles together

---
"    
  )
  
}



## -----------------------------------------------------------------------------
tidy_counts


## -----------------------------------------------------------------------------
counts_metadata


## -----------------------------------------------------------------------------
tidy_counts_meta <- inner_join(tidy_counts, counts_metadata, by = c("ENTREZ" = "ID"))


## -----------------------------------------------------------------------------
left_join(tidy_counts, counts_metadata, by = c("ENTREZ" = "ID"))


## -----------------------------------------------------------------------------
tidy_counts_expressed <- right_join(tidy_counts, counts_metadata, by = c("ENTREZ" = "ID"))
tidy_counts_expressed %>% tail()


## -----------------------------------------------------------------------------
semi_join(counts_metadata, tidy_counts, by = c("ID" = "ENTREZ") )


## -----------------------------------------------------------------------------

anti_join(counts_metadata, tidy_counts, by = c("ID" = "ENTREZ") )



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Quickly manipulate data with dplyr

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Quickly manipulate data with dplyr
---
"    
  )
  
}



## -----------------------------------------------------------------------------
select(tidy_counts_meta , counts)


## -----------------------------------------------------------------------------
select(tidy_counts_meta, counts, ENTREZ)


## -----------------------------------------------------------------------------
select(tidy_counts_meta, -Sample)


## -----------------------------------------------------------------------------
select(tidy_counts_meta, CellType:SYMBOL)


## -----------------------------------------------------------------------------
filter(tidy_counts_meta, Sample == 'CD34_1')


## -----------------------------------------------------------------------------
filter(tidy_counts_meta, Sample %in% c('CD34_1', 'ORTHO_1'))


## -----------------------------------------------------------------------------
filter(tidy_counts_meta, counts > 0)


## -----------------------------------------------------------------------------
arrange(tidy_counts_meta, counts)


## -----------------------------------------------------------------------------

arrange(tidy_counts_meta, CellType, desc(counts))


## -----------------------------------------------------------------------------
mutate(tidy_counts_meta, scale(counts))


## -----------------------------------------------------------------------------
mutate(tidy_counts_meta, count_zscore = scale(counts))


## -----------------------------------------------------------------------------
filter(tidy_counts_meta, counts == 0)


## -----------------------------------------------------------------------------
filter(tidy_counts_meta, counts == 0) %>%
  group_by(Sample) 


## -----------------------------------------------------------------------------
filter(tidy_counts_expressed, counts == 0) %>%
  group_by(Sample) %>%
  summarise(n())


## -----------------------------------------------------------------------------
tidy_counts_meta %>%
  group_by(ENTREZ) %>%
  summarise(counts_mean = mean(counts))


## -----------------------------------------------------------------------------
tidy_counts %>%
  group_by(ENTREZ, CellType) %>%
  summarise(counts_mean = mean(counts))


## -----------------------------------------------------------------------------
tidy_counts_meta %>%
  group_by(Sample) %>%
  filter(order(counts, decreasing=T) <= 3)


## -----------------------------------------------------------------------------
tidy_counts_meta %>% 
  filter(counts != 0) %>% 
  group_by(CellType, ENTREZ)


## -----------------------------------------------------------------------------
tidy_counts_meta %>% 
  filter(counts != 0) %>% 
  group_by(CellType, ENTREZ) %>%  
  filter(n()>1)


## -----------------------------------------------------------------------------
p <- tidy_counts_meta %>%
  group_by(ENTREZ, CellType) %>%
  summarise(counts_mean = mean(counts)) %>% 
  pivot_wider(names_from=CellType, values_from=counts_mean) %>%
  ggplot(aes(x=CD34, y=ORTHO)) + geom_point()



## -----------------------------------------------------------------------------
p 



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Outputting your tidy data

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Outputting your tidy data
---
"    
  )
  
}



## -----------------------------------------------------------------------------

write_delim(tidy_counts_meta, '../counts_with_metadata.csv', delim =',')

write_csv(tidy_counts_meta, '../counts_with_metadata.csv')


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Pattern matching with strings

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Pattern matching with strings
---
"    
  )
  
}



## -----------------------------------------------------------------------------

brc <- c("Tom", "Ji-Dung", "Matthew", "Wei", "Doug")


brc %>% str_sub(1, 3)


## -----------------------------------------------------------------------------
brc %>% str_sub(2, -2)


## -----------------------------------------------------------------------------

str_sub(brc, 2, -2) <- 'X'
brc


## -----------------------------------------------------------------------------
str_replace_all(brc, 'Matthew', 'Matt')
str_replace_all(brc, 'u', 'z' )


## -----------------------------------------------------------------------------
brc2 <- c("Tom  ", "  Ji  -Dung", "Matt   ", "Wei", "D o u g")

str_replace_all(brc2, ' ','' )


## -----------------------------------------------------------------------------
str_trim(brc2)



## -----------------------------------------------------------------------------

str_pad(brc2, width=10, side='left')


## -----------------------------------------------------------------------------

tidy_counts_meta %>% 
  pull(SYMBOL) %>%
  head()


tidy_counts_meta %>% 
  pull(SYMBOL) %>% 
  str_to_title() %>% 
  head()


## -----------------------------------------------------------------------------
tidy_counts_meta %>% 
  mutate(SYMBOL2 = str_to_title(SYMBOL))


## -----------------------------------------------------------------------------

tidy_counts_meta %>% 
  mutate(SYMBOL2 = str_to_title(SYMBOL)) %>% 
  mutate(SYMBOL3 = str_to_upper(SYMBOL2))


## -----------------------------------------------------------------------------

tidy_counts_meta %>% 
  pull(SYMBOL) %>% 
  str_detect('GAP') %>% head()

tidy_counts_meta %>% 
  pull(SYMBOL) %>% 
  str_detect('GAP') %>%
  filter(tidy_counts_meta, .)


## -----------------------------------------------------------------------------

tidy_counts_meta %>% 
  pull(SYMBOL) %>% 
  str_subset('GAP')



## -----------------------------------------------------------------------------

tidy_counts_meta %>% 
  pull(SYMBOL) %>% 
  str_count('GAP') %>% .[1:100]

tidy_counts_meta %>% 
  pull(SYMBOL) %>% 
  str_count('A') %>% .[1:100]

