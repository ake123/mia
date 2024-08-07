% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mia.R
\docType{data}
\name{peerj13075}
\alias{peerj13075}
\title{Skin microbial profiles 58 genetically unrelated individuals}
\format{
A TreeSummarizedExperiment with 674 features and 58 samples. The
rowData contains taxonomic information at kingdom, phylum, class, order,
family and genus level. The colData includes:

\describe{
\item{Sample}{sample ID}
\item{Geographical_location}{city where participant lives (Ahmednagar,
Pune and Nashik)}
\item{Gender}{participant's gender (Male or Female)}
\item{Age}{participant's age group (Middle_age, Adult and Elderly)}
\item{Diet}{participant's diet (Veg or Mixed)}
}
}
\usage{
data(peerj13075)
}
\description{
peerj13075 includes skin microbial profiles of 58 volunteers with multiple
factors. 16S r-RNA sequencing of V3-V4 regions was done to generate millions
of read using illumina platform. A standard bioinformatic and statistical
analysis done to explore skin bacterial diversity and its association with
age, diet, geographical locations. The authors investigated significant
association of skin microbiota with individual’s geographical location.
}
\references{
Potbhare, R., RaviKumar, A., Munukka, E., Lahti, L., & Ashma, R. (2022).
Skin microbiota diversity among genetically unrelated individuals of Indian origin.
PeerJ, 10, e13075. \url{https://doi.org/10.7717/peerj.13075}
Supplemental information includes OTU table and taxonomy table publicly-accessible from:
\url{https://www.doi.org/10.7717/peerj.13075/supp-1}
\url{https://www.doi.org/10.7717/peerj.13075/supp-2}
}
\seealso{
\code{\link{mia-datasets}}
}
\author{
Potbhare, R., et al.
}
\keyword{datasets}
