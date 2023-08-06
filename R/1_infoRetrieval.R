#' Extract PubMed or Europe PMC entry information for any query
#'
#' Combines functionality from two R-packages: "easyPubMed" and "europepmc" for reliable PubMed entry data extraction.
#' @name pubRetrieve
#' @param query Your typical PubMed query. Optimizing your query using the proper PubMed or Europe PMC syntax improves results.
#' @param size An integer representing the number of entries to retrieve. PubMed is limited to about 2000 articles, while Europe PMC can pull around 20,000.
#' @param database Define the PubMed database to retrieve from: "pubmed" or "pmc".
#' @keywords Information retrieval, PubMed, Europe PMC
#' @export
#' @examples
#' pubRetrieve(query = "Vape smoking AND toxicity", size = 2000, database = "pubmed")

pubRetrieve <- function(query, size, database) {

  #Dependencies:
  library(tidyverse)
  library(easyPubMed)
  library(europepmc)

      if(database == "pubmed") {

                  retrieved_info <- pubmed_retrieval(query = query, retmax = size)

      } else if(database == "pmc") {

                  retrieved_info <- europepmc_retrieval(query = query, retmax = size)

      } else {
            print("Error: database must equal 'pubmed' or 'pmc'")
      }#end if-'else if'-else loop

      return(retrieved_info)
}

#' PubMed multi-entry and retrieval function
#'
#' This function allows you to query and access key entry data from PubMed.
#' @name pubmed_retrieval
#' @param query Your typical PubMed query. Optimizing your query using the proper pubmed syntax improves results!
#' @param retmax Define how many entries you'd like to access for a query. Keep it between 500-1000 for fastest results. (Limit is around 2000)
#' @keywords Pubmed, PMC, Query, Retrieval
#' @examples
#' pubmed_retrieval(Query = "Vape smoking AND toxicity", retmax = 750)

pubmed_retrieval <- function(query, retmax) {

      #1 PMID Retrieval

      PMIDs <- get_pubmed_ids(query)
      print(paste(PMIDs$Count, "PMIDs retrieved. Fetching article information..."))

      #2 XML Format Record Download w/ PMID list
      articleInfo <- fetch_pubmed_data(PMIDs, retmax = retmax)

      #3 Convert XML to String
      #(this is to convert the content of each PubMed record to a character-class object)
      xmlToString <- articles_to_list(articleInfo)
      print(paste(length(xmlToString), "abstracts were retrieved. Creating output dataframe... (this may take a while)"))

      #4 Dataframe Retrieval
      stringToDF <- do.call(rbind,lapply(xmlToString, article_to_df, max_chars = -1, getAuthors = FALSE))
      return(stringToDF)
}

#' The Europe PMC search and retrieve function used for sourcing abstracts.
#'
#' This function allows you to query and access key entry data from European PubMed central.
#' @name europepmc_retrieval
#' @param query Your typical PubMed query. Optimizing your query using the proper pubmed syntax improves results!
#' @param retmax Define how many entries you'd like to access for a query. Keep it between 500-1000 for fastest results. (Limit is around 20,000)
#' @keywords Pubmed, EuropePMC, PMC, Query, Retrieval
#' @examples
#' europepmc_retrieval(Query = "Vape smoking AND toxicity", retmax = 750)

europepmc_retrieval <- function(query, retmax) {

            #1 PMID/PMC Retrieval
            europe_search <- epmc_search(query, output = 'raw', limit = retmax)

            #2 PMC or MEDLINE Retrieval
            retrieved_info <- keep(europe_search, function(x) (x[['source']] == 'MED') && !is.null(x[['abstractText']])) %>%
                  map_dfr(~ modify_at(.x, "journalInfo", flatten_dfr))

            #3 Standardizing column names
            names(retrieved_info)[names(retrieved_info) == "id"] <- "doc_id"
            names(retrieved_info)[names(retrieved_info) == "abstractText"] <- "abstract"

      return(retrieved_info)
}
