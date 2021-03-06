% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadSignatures.R
\name{loadCMSgenes}
\alias{loadCMSgenes}
\title{Load random-forest genes from CMS classifier,
stored in CMSclassifier::finalModel$importance, with their entrez identifiers}
\usage{
loadCMSgenes(geneSymbol = FALSE)
}
\arguments{
\item{geneSymbol, }{boolean: should entrez identifiers be converted to gene symbols; default is FALSE}
}
\value{
a matrix with CMS genes as rownames, and random forest importance scores per subtype;
        this is basically "CMSclassifier::finalModel$importance", with row names converted to
        gene symbols if geneSymbol=TRUE
}
\description{
Load random-forest genes from CMS classifier,
stored in CMSclassifier::finalModel$importance, with their entrez identifiers
}
\examples{
mat = loadCMSgenes(geneSymbol=FALSE)
print(head(mat))
}
\seealso{
\code{\link{loadSchlickerSignature}} and \code{\link{loadSadanandamSignature}}
}
