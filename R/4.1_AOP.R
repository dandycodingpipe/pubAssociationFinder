# KAF Rule AOP-wiki Fuzzy Matching
#' Compare KAF results to the AOP-wiki database to find notable key events, stressors, and biological processes present in your PubMed query.
#'
#' (In development) Uses locality-sensitive hashing to efficiently cluster and calculate the Jaccard similarity between mined rules and AOP-wiki values/classes. Values below 60% similarity are not considered matches.
#' @name matchAOP
#' @param rules Association rule dataframe that was mined using Abstract_ARM or easyKAF. Limit this parameter to 5000-20000 rules as deduplication is very computationally demanding.
#' @keywords Hashing, LSH, Locality sensitive hashing, deduplication, Jaccard Similarity, fuzzy matching, AOP
#' @export
#' @examples
#' KAFxAOP(rules)

KAFxAOP <- function(sample, AOPdownload) {

  results <- cbind(sample, AOPclass = NA, AOPvalue = NA)
  #and not FREE
if(length(sample$RHS) > 5000){
  lift_filtering <- results[order(-results$lift[1:5000]),]
} else {
  lift_filtering <- results[order(-results$lift),]
}

broken_down <- Rule_Concatenator(lift_filtering)

#Pre-downloaded aop-wiki quarterly backup data
doc <- AOPdownload

#XML Processing for AOP class/value extraction
doc <- AOP_XML_children_organizer(doc)

og_doc <- doc

#doc$value <- tolower(doc$value)

#Combine AOP and association-rule corpora
concat_df <- data.frame(corpora = c(doc$value, broken_down$toFuzzy))

# Locality Sensitive Hashing Package
library(textreuse)

#Defining hash table parameters
minhash <- minhash_generator(n = 6750, seed = 2553)

#Create hash codes for every object
Hashed <- TextReuseCorpus(text = concat_df$corpora, tokenizer = tokenize_ngrams, n = 1,
                          minhash_func = minhash, keep_tokens = TRUE,
                          progress = TRUE)

#Create buckets for storing potentially similar keys (hash codes)
Bucketing <- lsh(Hashed, bands = 450, progress = TRUE)

#Calculate and extract scores according to my threshold
final <- lsh_compare(lsh_candidates(Bucketing), Hashed, jaccard_similarity, progress = TRUE)
final = final[which(final$score >= 0.6),]


#The workflow of processing LSH results

#Removing extraneous characters and presenting results as coordinates
final$a = substr(final$a, 5, nchar(final$a))
final$b = substr(final$b, 5, nchar(final$b))

final$a <- as.numeric(final$a)
final$b <- as.numeric(final$b)

#Where are my AOP values
aoppies <- final[which(final$a <= 3956),]
aoppies$b = as.numeric(aoppies$b)

AOPcoordinate_list <- c()
Rulescoordinate_list <- c()
doc_val <- c()
aop_val <- c()
index_holder <- c()

momma <- data.frame()

for (k in 1:length(aoppies$a)) {
  if (aoppies$b[k] > 3956) {
    if (aoppies$score[k] >= 0.75) {
      index_holder <- c(index_holder, k)
      print(paste("PAIR", k, ":", concat_df[aoppies$a[k],], "+", concat_df[aoppies$b[k],], "SCORE:", aoppies$score[k]))
      print(which(doc$value == concat_df[aoppies$a[k],]))
      AOPcoordinate_list <- c(AOPcoordinate_list, which(doc$value == concat_df[aoppies$a[k],]))
      Rulescoordinate_list <- c(Rulescoordinate_list, which(broken_down$toFuzzy == concat_df[aoppies$b[k],]))

      # doc extract class info
      matching_indices <- which(doc$value == concat_df[aoppies$a[k],])
      doc_val <- c(doc_val, doc$AOP_Class[matching_indices])
      aop_val <- c(aop_val, doc$value[matching_indices])

      rule <- concat_df[aoppies$b[k],]
      whoa <- data.frame(rule = rule, aop_val = aop_val, doc_val = doc_val)
      momma <- rbind(momma, whoa)

      # Reset the vectors for the next iteration
      doc_val <- c()
      aop_val <- c()
      whoa <- data.frame()
    }
  }
}

length(unique(Rulescoordinate_list))
seerules <- lift_filtering[Rulescoordinate_list,]


original_index <- c()
for(i in 1:length(momma$rule)){
  original_index <- c(original_index, which(broken_down$toFuzzy == momma$rule[i]))
}

AOPresults <- momma
Rules <- seerules
return(list(AOPresults, Rules))
print("#####Done####")

}

# AOP-Wiki quarterly backup processing script
#' Converts AOP-wiki XML files into R data frames for further processing
#'
#' This function takes XML and converts
#' @name AOP_XML_children_organizer
#' @param doc XML file downloaded from AOP-wiki
#' @keywords Classification, AOP-wiki, AOP
#' @export
#' @examples
#' Word_Cleaner(data = association_rules_df)

AOP_XML_children_organizer <- function(doc) {

  library(XML)

  doc <- xmlRoot(xmlTreeParse(doc))

  class_frequencies <- table(names(doc))

  class_coordinates <- names(names(doc))

  unique_classes <- unique(class_coordinates)

  # Create a list containing the separated coordinates for every instance of every class

  my_list <- vector("list", length(unique_classes))
  names(my_list) <- unique_classes

  for(i in 1:length(unique_classes)){

    class_iteration_array <- which(class_coordinates == unique_classes[i])
    my_list[[i]] <- class_iteration_array

  }

  tmp = xmlSApply(doc, function(x) xmlSApply(x, xmlValue))


  tester <- data.frame(AOP_Class = 0, value = 0)

  for(h in 1:length(unique_classes)){

    label <- unique_classes[h]

    for(i in 1:length(my_list[[label]])){

      # Use switch case based on the value of the number variable
      result <- switch(label,
                       "chemical" = "preferred-name",
                       "biological-object" = "name",
                       "biological-process" = "name",
                       "biological-action" = "name",
                       "stressor" = "name",
                       "taxonomy" = "name",
                       "key-event" = "short-name",
                       "aop" = "short-name",
                       default = NULL)


      coordinate <- my_list[[label]][i]


      if (!is.null(result)) {
        text <- tmp[[coordinate]][result]

        new_row <- c(label, toString(text))
        tester <- rbind(tester, new_row)
      }

    }
  }
  return(tester)
}





