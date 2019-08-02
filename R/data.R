#' Nicotiana dataset
#'
#' A dataset containing a phylogenetic network and trait data for Nicotiana species
#'
#' The tree and data come from
#'
#' Chase M.W., Knapp S., Cox A.V., Clarkson J.J., Butsko Y., Joseph J., Savolainen V., and  Parokonny A.S. 2003. Molecular systematics, GISH and the origin of hybrid taxa in Nicotiana(Solanaceae). Annals of Botany 92: 107-127.
#'
#' Clarkson J.J., Lim K.Y., Kovarik A., Chase M.W., Knapp S. and Leitch A.R. 2005. Long-term genome diploidization I allopolyploid Nicotiana section Repandae(Solanaceae). New Phytologist 168:241-252.
#'
#' Komori T., Myers P.N., Yamada S., Kubo T., and Imaseki H. 2000. Comparative study of the Nicotiana species with respect to water deficit tolerance during early growth. Euphytica 116:121-130.
#'
#' @format A list with two items:
#' \describe{
#'   \item{phy.graph}{the phylogenetic network in ape::evonet format}
#'   \item{trait}{a vector of trait data}
#' }
"nicotiana"

#' Cichlid dataset
#'
#' A dataset containing a phylogenetic network and trait data for cichlid species
#'
#' The tree is made by doing a tree search with mitochondrial data from
#' Kobmuller, S., N. Duftner, K. M. Sefc, M. Aibara,M. Stipacek, M. Blanc, B. Egger, and C. Sturmbauer. 2007.Reticulate phylogeny of gastropod-shell-breeding cichlids from Lake Tanganyika: the result of repeated introgressive hybridization. BMC Evolutionary Biology 7:7.
#'
#' https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-7-7
#'
#' We then added hybridization events based on their cartoon Fig. 4: https://media.springernature.com/full/springer-static/image/art%3A10.1186%2F1471-2148-7-7/MediaObjects/12862_2006_Article_301_Fig4_HTML.jpg
#'
#' Hybridization events with solid lines (coeval events) were modeled as going from the later of the source or descendant nodes.
#'
#' Hybridization events with dotted lines, indicating ghost lineages, went from the MRCA of the source clade to the MRCA of the recipient taxon.
#'
#' Trait data comes from fishbase.
#'
#' @format A list with two items:
#' \describe{
#'   \item{phy.graph}{the phylogenetic network in ape::evonet format}
#'   \item{trait}{a vector of trait data}
#'   \item{final.se}{a vector of standard error}
#' }
"cichlid"
