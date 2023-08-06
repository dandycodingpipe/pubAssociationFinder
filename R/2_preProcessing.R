#' Reduce the size of your PubMed abstract texts
#'
#' Trim abstracts from the head and remove any entries without abstract text data.
#' @name abstractReduction
#' @param data PubMed entry data (see 'pubRetrieve')
#' @param fraction The decimal represenitng the fraction of each abstract to trim. We recommend (0.2) or 20%
#' @keywords Text reduction
#' @examples
#' abstractReduction(data = abstract_data, fraction = 0.2)

#1) reducing text size
      abstractReduction <- function(data, fraction){

            # Removes abstracts which contain "NA" instead of text
            NAs <- which(is.na(data$abstract))

                  if(is_empty(NAs)==TRUE) {
                        print('No corrupted abstracts in retrieval...')
                  } else {
                        data <- data[-NAs,]
                  }

            # Calculates the amount of characters to remove based on text size and removes a fraction
                  for(i in 1:length(data$abstract)) {

                        total <- nchar(data$abstract[i])
                        new_length <- total - round((fraction*total), digits = 0)
                        new_length <- as.integer(new_length)

                        data$abstract[i] <- str_trunc(data$abstract[i], new_length, 'left')
            }
            return(data)
      }#end abstract trimmer


      #' Parse retrieved abstract data using spaCy
      #'
      #' Parsed PubMed abstracts and extract potentially relevant information according to simple Parts-of-Speech or Language dependencies. For bio-entity extraction see '
      #' @name pubParse
      #' @param data PubMed entry data (see 'pubRetrieve')
      #' @param method Filtering method ('POS' or 'DEP'). Default: 'POS'
      #' @param composite 'Y' or defaults to null. Construct composite words using n-grams. This may help reduce trivial associations between words that form a composite word (e.g insulin -> resistance)
      #' @param venv_path REQUIRED: The path to your python venv. (see 'https://github.com/dandycodingpipe/KAFtool' for additional information)
      #' @param lang_model REQUIRED: The spaCy language model installed in your venv
      #' @param reduced_search The decimal represenitng the fraction of each abstract to trim. We recommend (0.2) or 20%
      #' @keywords Text parsing, sentence segmentation, tokenization, lemmatization, parts of speech, language dependencies, n-grams
      #' @export
      #' @examples
      #' pubParse(data = abstract_data, venv_path = "C:/Users/Chris/OneDrive/2023/Systox/venvJune19", lang_model = "en_core_web_lg", reduced_search = 0.2)

#2) one-and-done parser
      pubParse <- function(data, method, composite, venv_path, lang_model, reduced_search){
            library(tidyverse)
            library(spacyr)

            #1 Text reduction
            data <- abstractReduction(data, reduced_search)
            print(paste("1) Search reduction. Removing",reduced_search*100,"% from all abstract heads..."))

            #2 Establishing venv and R communication; sentence segmentation, tokenization, parts-of-speech, language dependcies, and lemmatization
            spacy_initialize(model = lang_model, virtualenv = venv_path)
            print("2) SpaCy initialized! Parsing the abstracts...")
            parsed_text <- spacy_parse(data$abstract, dependency = TRUE)

            #3) Removing sentences with negation modifiers
            parsed_text <- rm_negation(parsed_text)
            print('3) Sentences containing negation modifiers removed...')

            print('4) Filtering text...')
            # Stop word removal: parts-of-speech
            removePOS <- which(parsed_text$pos %in% c("PUNCT", "NUM", "SPACE"))
            parsed_text <- parsed_text[-removePOS,]

            if(is.null(method) | method == 'POS'){

              #4) POS filtering and lemma extraction
              POS <- which(parsed_text$pos %in% c("NOUN", "VERB", "ADJ"))
              parsed_text <- parsed_text[POS,]
              print(paste("filtered ",length(parsed_text$lemma),"potentially relevant lemma based on POS."))

              } else {

              #4) Dependency relationship filtering
              removeDPR <- which(parsed_text$dep_rel %in% c("det", "nummod", "punct", "aux", "prep", "cc", "ROOT", "auxpass",
                                                            "ccomp", "advmod", "xcomp", "relcl", "acl", "agent", "mark",
                                                            "expl", "pcomp", "intj", "meta", "csubj", "dative", "preconj",
                                                            "quantmod", "case", "predet", "csubjpass", "advcl", "nsubjpass",
                                                            "nsubj", "attr", "acomp", "poss", "oprd", "parataxis", "prt"))

              parsed_text <- parsed_text[-removeDPR,]
            }

            #5) Create composite words with processed text
            if(composite == "Y" | composite == "y"){

              library(RWeka)

              #Unparse text and regroup lemma for each document
              unparsed <- parsed_text %>%
                group_by(doc_id) %>%
                summarize(unparsed_text = paste(lemma, collapse = " ")) %>%
                ungroup()

              # Generate n-grams for each document that are 3 words large
              ngrams <- data.frame(ngrams = NGramTokenizer(unparsed, Weka_control(min = 1, max = 3)))

              # Detect duplicated n-grams
              duplicates <- data.frame(ngrams[duplicated(ngrams$ngrams), ])

              # Count the number of occurrences for each duplicated n-gram
              dup_counts <- as.data.frame(table(duplicates$ngrams))
              dup_counts <- dup_counts[which(dup_counts$Freq > 1),]

              # View the duplicated n-grams and their respective counts

              ####ngram integration and compound word preservation
              # Loop through each row of duplicate n-grams
              for (i in 1:nrow(dup_counts)) {
                ngram <- as.character(dup_counts$Var1[i])
                separator <- "_"  # or "-" or any other separator you prefer

                # Tokenize the n-gram into individual tokens
                tokens <- unlist(strsplit(ngram, " "))

                # Join the tokens with the separator in between
                separated_ngram <- paste(tokens, collapse = separator)

                # Replace instances of the duplicate n-gram with the separated version in the unparsed text
                unparsed$unparsed_text <- gsub(ngram, separated_ngram, unparsed$unparsed_text, fixed= TRUE)
              }

              reparsed <- spacy_parse(unparsed$unparsed_text)
              reparsed$lemma <- gsub("_", " ", reparsed$lemma)
              parsed_text <- reparsed
            }

            spacy_finalize()
            print('natural language processing complete!')

            return(parsed_text)
      }


      #' The function removes sentences possessing negations. Used in Text_Parser()
      #'
      #' This function removes negations in already-parsed abstract corpus.
      #' @name rm_negation
      #' @param Abstract_Parse Input your parsed abstracts. (typically obtained from the spacy_parse function)
      #' @keywords Text reduction, NLP, natural language processing, negation modifier
      #' @examples
      #' rm_negation(Abstract_Parse = parsed_abstracts_from_spacy_parse)

#3) negation removal function
      rm_negation <- function(Abstract_Parse){

            #3) Identifying Coordinates of Negative Modifiers

            Copy_Abstract_Parse <- Abstract_Parse
            Negation_Coord <- which(Copy_Abstract_Parse$dep_rel == "neg")
            print(paste("there are", length(Negation_Coord), "negations in the abstracts..."))

            #4) Function for Isolating Sentences Containing Negation Modifiers

            master_token_list <- c()
            print(paste("retrieving negation sentences..."))
            for(neg in 1:length(Negation_Coord)){

                  #identifying text_id and sentence_id corresponding to the negation coordinate
                  text_id <- Copy_Abstract_Parse$doc_id[Negation_Coord[neg]]

                  sentence_id <- Copy_Abstract_Parse$sentence_id[Negation_Coord[neg]]

                  #retrieving the coordinates for each token in the sentence containing negation
                  text_id_coordinates <- which(Copy_Abstract_Parse$doc_id == text_id)

                  sentence_id_coordinates <- which(Copy_Abstract_Parse$sentence_id[text_id_coordinates] == sentence_id)

                  token_ids <- text_id_coordinates[sentence_id_coordinates]

                  master_token_list <- c(master_token_list, token_ids)
            }
            print(paste("there are", length(master_token_list), "tokens associatied with negations in the abstracts. Removing..."))
            #Removing Sentences Containing Negations

            Copy_Abstract_Parse <- Copy_Abstract_Parse[-master_token_list,]
            return(Copy_Abstract_Parse)
      }
