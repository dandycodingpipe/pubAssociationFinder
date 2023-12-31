% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyKAF.R
\name{easyKAF}
\alias{easyKAF}
\title{Quickly mine PubMed literature from either PubMed or Europe PMC for signficant associations}
\usage{
easyKAF(query, database, venv_path, lang_model)
}
\arguments{
\item{query}{Your typical PubMed query. Optimizing your query using the proper PubMed or Europe PMC syntax improves results.}

\item{database}{Define the PubMed database to retrieve from: "pubmed" or "pmc". Currently, only "pubmed" articles can be used for bio-entity mining.}

\item{venv_path}{REQUIRED: The path to your python venv. (see 'https://github.com/dandycodingpipe/KAFtool' for additional information)}

\item{lang_model}{REQUIRED: The spaCy language model installed in your venv}
}
\description{
This method is for general literature exploration and limits users control for quick results. For more control over parameters, see 'pubRetrieve', 'pubParse', and 'pubMine'. For generating results that may be useful for translational research see 'bioKAF'. If this function does not work, you probably defined the venv wrong. You will need to restart R, and retry with the correct path.
}
\examples{
rules <- easyKAF(venv = "C:/Users/JohnDoe/venv/mar6", lang_model = "en_core_web_sm")
}
\keyword{ARM,}
\keyword{Association-rule}
\keyword{KAF,}
\keyword{KAFtool,}
\keyword{MeSH,}
\keyword{Systox}
\keyword{classification,}
\keyword{mining,}
\keyword{visualization,}
