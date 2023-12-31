


# MeSH Classification#

#This script makes exact word matches between the association rule consequents and MeSH's thesaurus



#test data set/ troubleshooting script
#test <- read.csv("sig_rules.csv")
#clean_words <- Word_Cleaner(test)
#exact <- MeSH_Mapper(test)

#word cleaning function
#' A crude consequent or RHS remover function for ambigious rules.
#'
#' This function removes the bracketing present in association-rules (e.g {consequent}) for better post-processing
#' @name Word_Cleaner
#' @param data Input your mined associations rules.
#' @keywords Post processing, Association-rule mining, ARM, classification
#' @examples
#' Word_Cleaner(data = association_rules_df)
library(tidyverse)
library(spacyr)

Word_Cleaner <- function(data) {

      #unique word list
      words <- unique(data$RHS)
      print(length(words))
      unique_rules <- as.character(words)

      #removes bracketing {RHS}
      for(i in 1:length(unique_rules)) {

            word <- toString(unique_rules[i])
            penul <- nchar(unique_rules[i])-1

            #substring the word by removing its first and last char
            unique_rules[i] <- substr(word,2,penul)

      }
      #removes problematic consequents that shouldnt exist anyways
      if(length(which(unique_rules=="%"))!=0) {
            unique_rules <- unique_rules[-which(unique_rules == "%")]
      }
      if(length(which(unique_rules=="="))!= 0) {
            unique_rules <- unique_rules[-which(unique_rules == "=")]
      }

      lilMatrix <- data.frame(unique_rules)
      return(lilMatrix)
}

#MeSH RDF API communicator for ancestral class extraction
#' A function that communicates with the Medical Subject Headings (MeSH) RDF API to classify consequenets (RHS) from mined association rules.
#'
#' This function is the classification workhouse which outputs the same association-rule dataframe you input, but with an additional column dedicated to the classification values.
#' @name MeSH_Mapper
#' @param word The data frame of rules you want to classify
#' @param removal A list of words you'd like to remove from the classification task (typically because they crash the process)
#' @keywords Post processing, Association-rule mining, ARM, classification
#' @examples
#' MeSH_Mapper(rules, removal = c("ml(-1","sub>2</sub","study","lead","±", "°"))
MeSH_Mapper <- function(word, removal){

removal = data.frame(unique_rules = removal)
words <- Word_Cleaner(word)
data <- words
coord_set <- c()
for(g in 1:length(removal$unique_rules)){
  coord <- which(words$unique_rules == removal[g,])
  coord_set = c(coord_set, coord)
}
if(!(is_empty(coord_set))){
words <- data.frame(unique_rules = words$unique_rules[-c(coord_set)])
}
      #this is a link that is used to interface with MeSH RDF API
      getsparq_left <- "https://id.nlm.nih.gov/mesh/sparql?query=PREFIX%20rdf%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23%3E%20PREFIX%20rdfs%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23%3E%20PREFIX%20xsd%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2001%2FXMLSchema%23%3E%20PREFIX%20owl%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2002%2F07%2Fowl%23%3E%20PREFIX%20meshv%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2Fvocab%23%3E%20PREFIX%20mesh%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F%3E%20PREFIX%20mesh2015%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2015%2F%3E%20PREFIX%20mesh2016%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2016%2F%3E%20PREFIX%20mesh2017%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2017%2F%3E%20%20%20SELECT%20%3Fd%20%3FdName%20%3Fc%20%3FcName%20%20FROM%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%3E%20%20%20%20WHERE%20%7B%20%20%20%20%20%20%3Fd%20a%20meshv%3ADescriptor%20.%20%20%3Fd%20meshv%3Aactive%201%20.%20%20%3Fd%20meshv%3Aconcept%20%3Fc%20.%20%20%3Fd%20rdfs%3Alabel%20%3FdName%20.%20%20%3Fc%20rdfs%3Alabel%20%3FcName%20%20FILTER(REGEX(%3FdName%2C%22"
      getsparq_right1 <- "%22%2C%22i%22)%20%7C%7C%20REGEX(%3FcName%2C%22"
      get_sqarq_right2 <-"%22%2C%22i%22))%20%20%20%20%7D%20%20ORDER%20BY%20%3Fd&format=JSON&inference=true&offset=0&limit=1000"

      SQL_token <- c()
      MeSH_Ancestry <- c()
      lilMatrix <- data.frame(words$unique_rules)

      #for every word,
      for(i in 1:length(words$unique_rules)) {
            print(paste0('Index ',i,": ",words$unique_rules[i]))
          # Replace spaces with '%20' if the string contains a space
            #words$unique_rules[i] <- ifelse(grepl(" ", words$unique_rules[i]), gsub(" ", "%20", words$unique_rules[i]), word$unique_rules[i])
            #words$unique_rules[i] <- ifelse(grepl("_", words$unique_rules[i]), gsub(" ", "%20", words$unique_rules[i]), word$unique_rules[i])
            #create the SQL query that returns the links for possible RDF id matches
            SQL_token_query_modifier <- paste(getsparq_left,words$unique_rules[i],getsparq_right1,words$unique_rules[i],get_sqarq_right2, sep = "")

            #get the possible word-match list and convert it to a usable data structure
            tokenID_Finder <- GET(SQL_token_query_modifier)
            tokenID_from_rdf <-fromJSON(rawToChar(tokenID_Finder$content))

            #capitalizing letters leaders to more frequent exact matches
            words$unique_rules[i] <- paste(toupper(substr(words$unique_rules[i], 1, 1)), substr(words$unique_rules[i], 2, nchar(words$unique_rules[i])), sep="")

            #print(paste("Rule:",unique_rules[i]))
            #check <- Match_Scoring(words$unique_rules[i], tokenID_from_rdf$results$bindings$dName$value,tokenID_from_rdf$results$bindings$cName$value)

            words$unique_rules[i]
            #check if any of the items are identical to the rule's string
            coord <- which(tokenID_from_rdf$results$bindings$cName$value == words$unique_rules[i])

            #this is a complicated if block that retrives a link if it exists (if i wanted to)
            if(length(coord) == 0){
                  coord <- which(tokenID_from_rdf$results$bindings$dName$value == words$unique_rules[i])


                  #secondary condition link if-statement
                  if(length(coord) == 0){
                              raw_link = 0
                  } else {
                        print(coord)
                        raw_link <- tokenID_from_rdf$results$bindings$d$value[coord]
                        }
            } else {
                  print(coord)
                  raw_link <- tokenID_from_rdf$results$bindings$d$value[coord]
            }

            raw_link <- toString(raw_link)
            #print(raw_link)

            if(nchar(raw_link) > 1){

                  RDF_ID <-substr(raw_link,28,nchar(raw_link))
                  #print(RDF_ID)
                  ancestry_L <- "https://id.nlm.nih.gov/mesh/sparql?query=PREFIX%20rdf%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23%3E%20PREFIX%20rdfs%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23%3E%20PREFIX%20xsd%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2001%2FXMLSchema%23%3E%20PREFIX%20owl%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2002%2F07%2Fowl%23%3E%20PREFIX%20meshv%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2Fvocab%23%3E%20PREFIX%20mesh%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F%3E%20PREFIX%20mesh2015%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2015%2F%3E%20PREFIX%20mesh2016%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2016%2F%3E%20PREFIX%20mesh2017%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2017%2F%3E%20%20SELECT%20%3FtreeNum%20%3FancestorTreeNum%20%3Fancestor%20%3Falabel%20FROM%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%3E%20%20WHERE%20%7B%20%20%20%20mesh%3A"
                  ancestry_R <- "%20meshv%3AtreeNumber%20%3FtreeNum%20.%20%20%20%20%3FtreeNum%20meshv%3AparentTreeNumber%2B%20%3FancestorTreeNum%20.%20%20%20%20%3Fancestor%20meshv%3AtreeNumber%20%3FancestorTreeNum%20.%20%20%20%20%3Fancestor%20rdfs%3Alabel%20%3Falabel%20%7D&format=JSON&inference=true&offset=0&limit=40"
                  ancestry_Query <- paste(ancestry_L,RDF_ID,ancestry_R, sep = '')
                  ancestry_Finder <- GET(ancestry_Query)

                  data <- fromJSON(rawToChar(ancestry_Finder$content))

                  d1 <-substr(data$results$bindings$ancestorTreeNum$value[1],28,nchar(data$results$bindings$ancestorTreeNum$value[1]))

                  d1 <- substr(d1,1,1)
                  if(is_empty(d1)== TRUE){
                        d1 <- "NA"
                  }
                  print(paste("Rule:",words$unique_rules[i], " | Class:",d1))
                  MeSH_Ancestry <- c(MeSH_Ancestry,d1)
            } else {
                        MeSH_Ancestry <- c(MeSH_Ancestry,"NA")
            }

            #print(MeSH_Ancestry)
      }

      lilMatrix <- cbind(words$unique_rules,MeSH_Ancestry)
      lilMatrix <- data.frame(lilMatrix)
      lilMatrix <- MeSH_filter(lilMatrix)


      return(lilMatrix)
}

#Recreate the bracketing in association rules
#' This function simply rebuilds the association rules that the MeSH cleaner had removed.
#'
#' Used for MeSH_Cleaner.
#' @name MeSH_filter
#' @param MeSH_data The data frame of rules you want to re-structure
#' @keywords Post processing, Association-rule mining, ARM, classification, MeSH
#' @examples
#' MeSH_filter(MeSH_data)
MeSH_filter <- function(MeSH_data){

      for( i in 1:length(MeSH_data$V1)){
            MeSH_data$V1[i] <- tolower(MeSH_data$V1[i])
            MeSH_data$V1[i] <- paste("{", MeSH_data$V1[i] ,"}",sep = '')
      }
      return(MeSH_data)
}

#final output
#' The KAF classification work-horse
#'
#' This function handles the classification suite of functions and outputs cleaned and classified association rules
#' @name matchMeSH
#' @param raw_rules The data frame of rules you want to classify
#' @param removal A list of words that crash the classifier. These words are removed and not classified.
#' @keywords Post processing, Association-rule mining, ARM, classification, MeSH
#' @export
#' @examples
#' matchMeSH(rules, removal = c("ml(-1","sub>2</sub","study","lead","±", "°", "confidence", "-", "%", ">", "sub>50</sub","müllerian", "#"))

matchMeSH <- function(raw_rules, removal){

    library(httr)
    library(jsonlite)
    library(dplyr)

      classifications <- MeSH_Mapper(raw_rules,removal)

      colnames(classifications)[1] ="RHS"

      final_Df <- left_join(raw_rules, classifications)

      return(final_Df)
}


#' Clean association rules
#'
#' Before fuzzy matching we have to pre-process the rules by concatenating them into a single string
#' @name Rule_Concatenator
#' @param rules The data frame of rules you want to concatenate
#' @keywords Post processing, Association-rule mining, ARM, classification, MeSH
#' @export
#' @examples Rule_Concatenator(rules)
Rule_Concatenator <- function(rules) {
  LHS <- as.character(rules$LHS)
  RHS <- as.character(rules$RHS)

  LHSpenul <- nchar(LHS) - 1
  RHSpenul <- nchar(RHS) - 1

  toFuzzy <- paste(substr(LHS, 2, LHSpenul), substr(RHS, 2, RHSpenul), sep = " ")

  return(data.frame(toFuzzy, rules$RHS, stringsAsFactors = FALSE))
}

