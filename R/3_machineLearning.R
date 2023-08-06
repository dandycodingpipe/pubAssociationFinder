

#' Mine association-rules on parsed PubMed abstract data
#'
#' The output is an R dataframe containing association rules and associated metrics. Outputs can vary greatly depending on the parsing method.
#' @name pubMine
#' @param data Parsed data (See pubParse)
#' @param min_supp Define the minimum support threshold for the rules (e.g 0.01 or 0.10)
#' @param min_conf Define the minimum confidence threshold for the rules (e.g 0.50 or 0.90)
#' @param min_p Define the minimum p-value threshold for the rules (e.g 0.05 or 0.0.005)
#' @keywords Processing algorithm, association-rules mining, statistics, arules,
#' @export
#' @examples
#' pubMine(data <- data_from_Text_Parser, min_supp = 0.01, min_conf = 0.75, min_p = 0.005)

pubMine <- function(data, min_supp, min_conf, min_p ){
  library(arules)

  print("Converting pasrsed abstracts to item/transaction format...")
      data$doc_id <- as.factor(data$doc_id)
      data$lemma <- as.factor(data$lemma)
      txns <- as(split(data$lemma,data$doc_id), "transactions")

  print("Initiating apriori algorithm...")
  rules <- apriori(txns, parameter = list(supp = min_supp, conf = min_conf, target = 'rules'), control = list(memopt = TRUE))

  #Statistical filtering
  print("Removing rules that do not meet p-value thresholds and type-1 errors...")
  table <- is.significant(rules, txns, method = "Fisher", alpha = min_p, adjust = "fdr")
  falses <- which(table == FALSE)
  rules <- rules[-falses,]

  print("Formatting rules into R-data frame and removing special symbols of n:...")
  df_rules <- DATAFRAME(rules)
    if(length(which(df_rules$RHS=="{-}"))!= 0) {
      print(length(which(df_rules$RHS=="{-}")))
      df_rules <- df_rules[-which(df_rules$RHS == "{-}"),]
    }

  if(length(which(df_rules$RHS=="{%}"))!= 0) {
    print(length(which(df_rules$RHS=="{%}")))
      df_rules <- df_rules[-which(df_rules$RHS== "{%}"),]
  }

  if(length(which(df_rules$RHS=="{=}"))!= 0) {
    print(length(which(df_rules$RHS=="{=}")))
    df_rules <- df_rules[-which(df_rules$RHS== "{=}"),]
  }
  print("Done!")
  return(df_rules)
}





