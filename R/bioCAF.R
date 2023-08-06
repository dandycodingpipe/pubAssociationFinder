# PubTator Interface
#' Mine bioconcept associations by combining PubTator and KAF
#'
#' Bioconcepts associations this function can produce are chemical-gene, gene-gene, or disease-chemical/disease-gene.
#' @name bioCAF
#' @param data Generated using PubTator
#' @param choice Control which associations to mine and retrieve for analysis: "All", Chemical", "Disease", "Gene", "Systox" (Chemical and disease associations)
#' @param support To control computational load, support is modifiable according to each query and user needs. If you are generating too many results for your liking, you can increase this support threshold which represents the percentage of (=<2000) abstracts
#' @param venv_path REQUIRED: The path to your python venv. (see 'https://github.com/dandycodingpipe/KAFtool' for additional information)
#' @param lang_model REQUIRED: The spaCy language model installed in your venv
#' @keywords KAF, PubTator, PubMed
#' @export
#' @examples
#' bioCAF(data, "Disease", 0.001, "C:/Users/Chris/OneDrive/2023/Systox/venvJune19", "en_core_web_lg")

bioCAF <- function(data, choice, support, venv_path, lang_model) {
  entryData <- data[[1]]
  BioAnnotations <- data[[2]]

  for (j in 1:length(entryData$pmid)) {

    data_text <- BioAnnotations[j]

    library(stringr)


    # Read the data into a data frame
    df <- read.delim(text = data_text, header = FALSE, stringsAsFactors = FALSE)
    og_length <- length(df$V6)
    if(og_length == 0){
     next
    }
    df <- df[-c(1,2),]
    print(j)
    #print(paste((which(duplicated(df$V6) == FALSE)/og_length)*100, "% reduced"))
    df <- df[which(duplicated(df$V6) == FALSE),]

    # Filter out the bioconcepts
    Chemwords <- unique(df$V4[which(df$V5 == "Chemical")])
    if(!is_empty(Chemwords)){
    Chemwords <- paste0("c_", Chemwords)
    }
    Geneword <- unique(df$V4[which(df$V5 == "Gene")])
    if(!is_empty(Geneword)){
    Geneword <- paste0("g_", Geneword)
    }
    Mutation <- unique(df$V4[which(df$V5 == "Mutation")])
    if(!is_empty(Mutation)){
    Mutation <- paste0("m_", Mutation)
    }
    Diseasewords <- unique(df$V4[which(df$V5 == "Disease")])
    if(!is_empty(Diseasewords)){
      Diseasewords <- paste0("d_", Diseasewords)
    }


    if(choice == "All" | choice == "all"){
      bioconcepts <- c(Chemwords, Geneword, Diseasewords, Mutation)

    # Disease-gene associations
    }else if(choice == "Disease" | choice == "disease"){
      bioconcepts <- c(Diseasewords, Geneword, Mutation)
      print("here")
    }
    # Chemical-gene associations
    else if(choice == "Chemical" | choice == "chemical"){
      bioconcepts <- c(Chemwords, Geneword, Mutation)
    }
    # Gene-gene associations
    else if(choice == "Gene" | choice == "gene"){
      bioconcepts <- c(Geneword)
      print('here')
    }
    # Chemical-disease associations
    else if(choice == "Systox" | choice == "systox"){
      bioconcepts <- c(Chemwords, Diseasewords)
      print('here')
    }


     bioconcepts <- add_underscores_to_compound_words(bioconcepts)

    # Concatenate bioconcepts as a single string and store it in the Total_Data_Copy data frame
    entryData$bioannotations[j] <- paste(tolower(bioconcepts), collapse = " ")
  }

  #NLP / Machine Learning Quick-Suite

  library(spacyr)
  spacy_initialize(model = lang_model, virtualenv = venv)
  parsed <- spacy_parse(entryData$bioannotations)
  cleaned <- parsed[-c(which(parsed$token <= 1)),]

  library(arules)
  V2associations <- pubMine(cleaned, support, 0.6, 0.05)
  return(V2associations)
  }

#' Replace bioentity noun phrase spacing with underscores to preserve structure
#'
#' This function is not intended for users
#' @name add_underscores_to_compound_words
#' @param character_list List of named entities in character type
  add_underscores_to_compound_words <- function(character_list) {
    # Use gsub to replace spaces with underscores in each element of the character list
    updated_list <- gsub(" ", "_", character_list)
    updated_list <- gsub("-", "_", updated_list)


    return(updated_list)
}




