#' Motif Database
#'
#' A compilation of DNA motifs from various sources.
#'
#' @format A data frame with 270340 rows and 9 columns:
#' \describe{
#'   \item{motif}{Character. The name or identifier of the motif.}
#'   \item{pos}{Integer. The position within the motif.}
#'   \item{A}{Numeric. The score for adenine at this position.}
#'   \item{C}{Numeric. The score for cytosine at this position.}
#'   \item{G}{Numeric. The score for guanine at this position.}
#'   \item{T}{Numeric. The score for thymine at this position.}
#'   \item{dataset}{Character. The source dataset of the motif.}
#'   \item{motif_orig}{Character. The name of the motif without the database prefix.}
#'   \item{cluster}{Character. For SCENIC motifs, the cluster to which the motif belongs.}
#' }
#'
#' @references
#' \describe{
#' \item{HOMER: }{Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432}
#' \item{JASPAR: }{Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I, Berhanu Lemma R, Turchi L, Blanc-Mathieu R, Lucas J, Boddie P, Khan A, Manosalva Pérez N, Fornes O, Leung TY, Aguirre A, Hammal F, Schmelter D, Baranasic D, Ballester B, Sandelin A, Lenhard B, Vandepoele K, Wasserman WW, Parcy F, and Mathelier A JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles Nucleic Acids Res. 2022 Jan 7;50(D1):D165-D173.; doi: 10.1093/nar/gkab1113}
#' \item{JOLMA: }{Jolma, A., Yin, Y., Nitta, K. et al. DNA-dependent formation of transcription factor pairs alters their binding specificity. Nature 534, S15–S16 (2016). \url{https://doi.org/10.1038/nature18912}}
#' \item{HOCOMOCO: }{Ivan V. Kulakovskiy; Ilya E. Vorontsov; Ivan S. Yevshin; Ruslan N. Sharipov; Alla D. Fedorova; Eugene I. Rumynskiy; Yulia A. Medvedeva; Arturo Magana-Mora; Vladimir B. Bajic; Dmitry A. Papatsenko; Fedor A. Kolpakov; Vsevolod J. Makeev: HOCOMOCO: towards a complete collection of transcription factor binding models for human and mouse via large-scale ChIP-Seq analysis. Nucl. Acids Res., Database issue, gkx1106 (11 November 2017). \url{https://doi.org/10.1093/nar/gkx1106}}
#' \item{SCENIC: }{Aibar, S., González-Blas, C., Moerman, T. et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods 14, 1083–1086 (2017). \url{https://doi.org/10.1038/nmeth.4463}}
#' \item{SCENIC+: }{Bravo González-Blas, C., De Winter, S., Hulselmans, G. et al. SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. Nat Methods 20, 1355–1367 (2023). \url{https://doi.org/10.1038/s41592-023-01938-4}}
#' }
#'
#' @examples
#' data(motif_db)
#' head(motif_db)
"motif_db"
