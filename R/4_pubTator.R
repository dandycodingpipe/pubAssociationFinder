#Extract bioconcepts from PubMed abstracts
#' Retrieve bio-entities in PubMed data using PubTator
#'
#' This function uses PubTator-annotated articles to retrieve pre-annotated articles. These results should be passed onto bioKAF for further processing (See bioKAF).
#' @name PubTator
#' @param query Input a PubMed query (see https://pubmed.ncbi.nlm.nih.gov/help/#how-do-i-search-pubmed)
#' @keywords KAF, PubTator, PubMed
#' @export
#' @examples
#' PubTator("dioxin toxicity","C:/Users/Chris/OneDrive/2023/Systox/venvJune19", "en_core_web_lg" )
PubTator <- function(query) {

  library(easyPubMed)
  library(bioTM)
  library(stringr)
  library(dplyr)
  library(httr)
  library(tidyverse)

  Total_Data <- info_retrieval(query, 2000, "pubmed")

  BioAnnotations <- c()
  REDUCTION <- c()

  for (i in 1:length(Total_Data$pmid)) {
    id <- as.character(Total_Data$pmid[i])
    r_url_START <- "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/BioConcept/"
    r_url_END <- "/PubTator/"
    r_url <- paste0(r_url_START, id, r_url_END)

    response <- GET(r_url)

    result <- rawToChar(response$content)

    if (result == "") {
      toss <- which(Total_Data$pmid == Total_Data$pmid[i])
      REDUCTION <- c(REDUCTION, toss)
    } else {
      BioAnnotations <- c(BioAnnotations, result)
      print(paste("Document:",i))
    }
  }

  Total_Data_Copy <- Total_Data[-REDUCTION,]
  return(list(Total_Data_Copy, BioAnnotations))
}
