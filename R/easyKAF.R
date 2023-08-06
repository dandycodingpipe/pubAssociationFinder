#' Quickly mine PubMed literature from either PubMed or Europe PMC for signficant associations
#'
#' This method is for general literature exploration and limits users control for quick results. For more control over parameters, see 'pubRetrieve', 'pubParse', and 'pubMine'. For generating results that may be useful for translational research see 'bioKAF'. If this function does not work, you probably defined the venv wrong. You will need to restart R, and retry with the correct path.
#' @name easyKAF
#' @param query Your typical PubMed query. Optimizing your query using the proper PubMed or Europe PMC syntax improves results.
#' @param database Define the PubMed database to retrieve from: "pubmed" or "pmc". Currently, only "pubmed" articles can be used for bio-entity mining.
#' @param venv_path REQUIRED: The path to your python venv. (see 'https://github.com/dandycodingpipe/KAFtool' for additional information)
#' @param lang_model REQUIRED: The spaCy language model installed in your venv
#' @keywords KAF, Association-rule mining, ARM, classification, MeSH, visualization, KAFtool, Systox
#' @export
#' @examples
#' rules <- easyKAF(venv = "C:/Users/JohnDoe/venv/mar6", lang_model = "en_core_web_sm")
easyKAF <- function(query, database, venv_path, lang_model){

      # 1. Retrieval
        retrieved <- pubRetrieve(query, 1500, database)

      # 2. Natural language processing
        preProcessed <- pubParse(retrieved, method = "POS", composite = "n", venv_path, lang_model, 0.2)

      # 3. Unsupervised machine learning
        rules <- pubMine(preProcessed, 0.01, 0.75, 0.005)

      # 3.1 Pruning large-rule set
        filt_rules <- which(rules$lift <= 2)
        rules <- rules[-filt_rules,]

        return(rules)
}


