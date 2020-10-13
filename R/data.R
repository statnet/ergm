#  File R/data.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' Two versions of an E. Coli network dataset
#' 
#' This network data set comprises two versions of a biological network in
#' which the nodes are operons in \emph{Escherichia Coli} and a directed edge
#' from one node to another indicates that the first encodes the transcription
#' factor that regulates the second.
#' 
#' The network object \code{ecoli1} is directed, with 423 nodes and 519 arcs.
#' The object \code{ecoli2} is an undirected version of the same network, in
#' which all arcs are treated as edges and the five isolated nodes (which
#' exhibit only self-regulation in \code{ecoli1}) are removed, leaving 418
#' nodes.
#' 
#' @usage
#' data(ecoli)
#' @docType data
#' @name ecoli
#' @aliases ecoli ecoli1 ecoli2
#' @section Licenses and Citation: When publishing results obtained using this
#' data set, the original authors (Salgado et al, 2001; Shen-Orr et al, 2002)
#' should be cited, along with this \code{R} package.
#' @references
#' 
#' Salgado et al (2001), Regulondb (version 3.2): Transcriptional Regulation
#' and Operon Organization in Escherichia Coli K-12, \emph{Nucleic Acids
#' Research}, 29(1): 72-74.
#' 
#' Shen-Orr et al (2002), Network Motifs in the Transcriptional Regulation
#' Network of Escerichia Coli, \emph{Nature Genetics}, 31(1): 64-68.
#' 
#' %Saul and Filkov (2007)
#' 
#' %Hummel et al (2010)
#' @source The data set is based on the RegulonDB network (Salgado et al, 2001)
#' and was modified by Shen-Orr et al (2002).
#' @keywords datasets
NULL

#' Faux desert High School as a network object
#' 
#' This data set represents a simulation of a directed in-school friendship
#' network.  The network is named faux.desert.high.
#' 
#' 
#' @usage
#' data(faux.desert.high)
#' @docType data
#' @name faux.desert.high
#' @format \code{faux.desert.high} is a \code{\link[network]{network}} object
#' with 107 vertices (students, in this case) and 439 directed edges
#' (friendship nominations). To obtain additional summary information about it,
#' type \code{summary(faux.desert.high)}.
#' 
#' The vertex attributes are \code{Grade}, \code{Sex}, and \code{Race}. The
#' \code{Grade} attribute has values 7 through 12, indicating each student's
#' grade in school.  The \code{Race} attribute is based on the answers to two
#' questions, one on Hispanic identity and one on race, and takes six possible
#' values: White (non-Hisp.), Black (non-Hisp.), Hispanic, Asian (non-Hisp.),
#' Native American, and Other (non-Hisp.)
#' @section Licenses and Citation: If the source of the data set does not
#' specified otherwise, this data set is protected by the Creative Commons
#' License \url{https://creativecommons.org/licenses/by-nc-nd/2.5/}.
#' 
#' When publishing results obtained using this data set, the original authors
#' (Resnick et al, 1997) should be cited. In addition this package should be
#' cited as:
#' 
#' Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau, and
#' Martina Morris. 2003 \emph{statnet: Software tools for the Statistical
#' Modeling of Network Data} \cr \url{https://statnet.org}.
#' @seealso \code{\link[network]{network}},
#' \code{\link[network]{plot.network}}, \code{\link{ergm}},
#' \code{\link{faux.desert.high}}, \code{\link{faux.mesa.high}},
#' \code{\link{faux.magnolia.high}}
#' @references
#' 
#' Resnick M.D., Bearman, P.S., Blum R.W. et al. (1997). \emph{Protecting
#' adolescents from harm. Findings from the National Longitudinal Study on
#' Adolescent Health}, \emph{Journal of the American Medical Association}, 278:
#' 823-32.
#' @source The data set is simulation based upon an ergm model fit to data from
#' one school community from the AddHealth Study, Wave I (Resnick et al.,
#' 1997). It was constructed as follows:
#' 
#' The school in question (a single school with 7th through 12th grades) was
#' selected from the Add Health "structure files."  Documentation on these
#' files can be found here:
#' \url{https://addhealth.cpc.unc.edu/documentation/codebooks/}.
#' 
#' The stucture file contains directed out-ties representing each instance of a
#' student who named another student as a friend.  Students could nominate up
#' to 5 male and 5 female friends. Note that registered students who did not
#' take the AddHealth survey or who were not listed by name on the schools'
#' student roster are not included in the stucture files.  In addition, we
#' removed any students with missing values for race, grade or sex.
#' 
#' The following \code{\link{ergm}} model was fit to the original data:
#' 
#' \preformatted{ desert.fit <- ergm(original.net ~ edges + mutual +
#' absdiff("grade") + nodefactor("race", base=5) + nodefactor("grade", base=3)
#' + nodefactor("sex") + nodematch("race", diff = TRUE) + nodematch("grade",
#' diff = TRUE) + nodematch("sex", diff = FALSE) + idegree(0:1) + odegree(0:1)
#' + gwesp(0.1,fixed=T), constraints = ~bd(maxout=10), control =
#' control.ergm(MCMLE.steplength = .25, MCMC.burnin = 100000, MCMC.interval =
#' 10000, MCMC.samplesize = 2500, MCMLE.maxit = 100), verbose=T) }
#' 
#' Then the faux.desert.high dataset was created by simulating a single network
#' from the above model fit:
#' 
#' \preformatted{ faux.desert.high <- simulate(desert.fit, nsim=1, burnin=1e+8,
#' constraint = "edges") }
#' @keywords datasets
NULL

#' Faux dixon High School as a network object
#' 
#' This data set represents a simulation of a directed in-school friendship
#' network.  The network is named faux.dixon.high.
#' 
#' 
#' @usage
#' data(faux.dixon.high)
#' @docType data
#' @name faux.dixon.high
#' @format \code{faux.dixon.high} is a \code{\link[network]{network}} object
#' with 248 vertices (students, in this case) and 1197 directed edges
#' (friendship nominations). To obtain additional summary information about it,
#' type \code{summary(faux.dixon.high)}.
#' 
#' The vertex attributes are \code{Grade}, \code{Sex}, and \code{Race}. The
#' \code{Grade} attribute has values 7 through 12, indicating each student's
#' grade in school.  The \code{Race} attribute is based on the answers to two
#' questions, one on Hispanic identity and one on race, and takes six possible
#' values: White (non-Hisp.), Black (non-Hisp.), Hispanic, Asian (non-Hisp.),
#' Native American, and Other (non-Hisp.)
#' @section Licenses and Citation: If the source of the data set does not
#' specified otherwise, this data set is protected by the Creative Commons
#' License \url{https://creativecommons.org/licenses/by-nc-nd/2.5/}.
#' 
#' When publishing results obtained using this data set, the original authors
#' (Resnick et al, 1997) should be cited. In addition this package should be
#' cited as:
#' 
#' Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau, and
#' Martina Morris. 2003 \emph{statnet: Software tools for the Statistical
#' Modeling of Network Data} \cr \url{https://statnet.org}.
#' @seealso \code{\link[network]{network}},
#' \code{\link[network]{plot.network}}, \code{\link{ergm}},
#' \code{\link{faux.desert.high}}, \code{\link{faux.mesa.high}},
#' \code{\link{faux.magnolia.high}}
#' @references
#' 
#' Resnick M.D., Bearman, P.S., Blum R.W. et al. (1997). \emph{Protecting
#' adolescents from harm. Findings from the National Longitudinal Study on
#' Adolescent Health}, \emph{Journal of the American Medical Association}, 278:
#' 823-32.
#' @source The data set is simulation based upon an ergm model fit to data from
#' one school community from the AddHealth Study, Wave I (Resnick et al.,
#' 1997). It was constructed as follows:
#' 
#' The school in question (a single school with 7th through 12th grades) was
#' selected from the Add Health "structure files."  Documentation on these
#' files can be found here:
#' \url{https://addhealth.cpc.unc.edu/documentation/codebooks/}.
#' 
#' The stucture file contains directed out-ties representing each instance of a
#' student who named another student as a friend.  Students could nominate up
#' to 5 male and 5 female friends. Note that registered students who did not
#' take the AddHealth survey or who were not listed by name on the schools'
#' student roster are not included in the stucture files.  In addition, we
#' removed any students with missing values for race, grade or sex.
#' 
#' The following \code{\link{ergm}} model was fit to the original data:
#' 
#' \preformatted{ dixon.fit <- ergm(original.net ~ edges + mutual +
#' absdiff("grade") + nodefactor("race", base=5) + nodefactor("grade", base=3)
#' + nodefactor("sex") + nodematch("race", diff = TRUE) + nodematch("grade",
#' diff = TRUE) + nodematch("sex", diff = FALSE) + idegree(0:1) + odegree(0:1)
#' + gwesp(0.1,fixed=T), constraints = ~bd(maxout=10), control =
#' control.ergm(MCMLE.steplength = .25, MCMC.burnin = 100000, MCMC.interval =
#' 10000, MCMC.samplesize = 2500, MCMLE.maxit = 100), verbose=T) }
#' 
#' Then the faux.dixon.high dataset was created by simulating a single network
#' from the above model fit:
#' 
#' \preformatted{ faux.dixon.high <- simulate(dixon.fit, nsim=1, burnin=1e+8,
#' constraint = "edges") }
#' @keywords datasets
NULL

#' Goodreau's Faux Magnolia High School as a network object
#' 
#' This data set represents a simulation of an in-school friendship network.
#' The network is named faux.magnolia.high because the school commnunities on
#' which it is based are large and located in the southern US.
#' 
#' 
#' @usage
#' data(faux.magnolia.high)
#' @docType data
#' @name faux.magnolia.high
#' @format \code{faux.magnolia.high} is a \code{\link[network]{network}} object
#' with 1461 vertices (students, in this case) and 974 undirected edges (mutual
#' friendships). To obtain additional summary information about it, type
#' \code{summary(faux.magnolia.high)}.
#' 
#' The vertex attributes are \code{Grade}, \code{Sex}, and \code{Race}. The
#' \code{Grade} attribute has values 7 through 12, indicating each student's
#' grade in school.  The \code{Race} attribute is based on the answers to two
#' questions, one on Hispanic identity and one on race, and takes six possible
#' values: White (non-Hisp.), Black (non-Hisp.), Hispanic, Asian (non-Hisp.),
#' Native American, and Other (non-Hisp.)
#' @section Licenses and Citation: If the source of the data set does not
#' specified otherwise, this data set is protected by the Creative Commons
#' License \url{https://creativecommons.org/licenses/by-nc-nd/2.5/}.
#' 
#' When publishing results obtained using this data set, the original authors
#' (Resnick et al, 1997) should be cited. In addition this package should be
#' cited as:
#' 
#' Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau, and
#' Martina Morris. 2003 \emph{statnet: Software tools for the Statistical
#' Modeling of Network Data} \cr \url{https://statnet.org}.
#' @seealso \code{\link[network]{network}},
#' \code{\link[network]{plot.network}}, \code{\link{ergm}},
#' \code{\link{faux.mesa.high}}
#' @references
#' 
#' Resnick M.D., Bearman, P.S., Blum R.W. et al. (1997). \emph{Protecting
#' adolescents from harm. Findings from the National Longitudinal Study on
#' Adolescent Health}, \emph{Journal of the American Medical Association}, 278:
#' 823-32.
#' @source The data set is based upon a model fit to data from two school
#' communities from the AddHealth Study, Wave I (Resnick et al., 1997). It was
#' constructed as follows:
#' 
#' The two schools in question (a junior and senior high school in the same
#' community) were combined into a single network dataset.  Students who did
#' not take the AddHealth survey or who were not listed on the schools' student
#' rosters were eliminated, then an undirected link was established between any
#' two individuals who both named each other as a friend.  All missing race,
#' grade, and sex values were replaced by a random draw with weights determined
#' by the size of the attribute classes in the school.
#' 
#' The following \code{\link{ergm}} model was fit to the original data:
#' 
#' \preformatted{ magnolia.fit <- ergm (magnolia ~ edges +
#' nodematch("Grade",diff=T) + nodematch("Race",diff=T) +
#' nodematch("Sex",diff=F) + absdiff("Grade") + gwesp(0.25,fixed=T),
#' burnin=10000, interval=1000, MCMCsamplesize=2500, maxit=25,
#' control=control.ergm(steplength=0.25)) }
#' 
#' Then the faux.magnolia.high dataset was created by simulating a single
#' network from the above model fit:
#' 
#' \preformatted{ faux.magnolia.high <- simulate (magnolia.fit, nsim=1,
#' burnin=100000000, constraint = "edges") }
#' @keywords datasets
NULL

#' Goodreau's Faux Mesa High School as a network object
#' 
#' This data set (formerly called \dQuote{fauxhigh}) represents a simulation of
#' an in-school friendship network.  The network is named \code{faux.mesa.high}
#' because the school commnunity on which it is based is in the rural western
#' US, with a student body that is largely Hispanic and Native American.
#' 
#' 
#' @usage
#' data(faux.mesa.high)
#' @docType data
#' @name faux.mesa.high
#' @aliases faux.mesa.high fauxhigh
#' @format \code{faux.mesa.high} is a \code{\link[network]{network}} object
#' with 205 vertices (students, in this case) and 203 undirected edges (mutual
#' friendships).  To obtain additional summary information about it, type
#' \code{summary(faux.mesa.high)}.
#' 
#' The vertex attributes are \code{Grade}, \code{Sex}, and \code{Race}. The
#' \code{Grade} attribute has values 7 through 12, indicating each student's
#' grade in school.  The \code{Race} attribute is based on the answers to two
#' questions, one on Hispanic identity and one on race, and takes six possible
#' values: White (non-Hisp.), Black (non-Hisp.), Hispanic, Asian (non-Hisp.),
#' Native American, and Other (non-Hisp.)
#' @section Licenses and Citation: If the source of the data set does not
#' specified otherwise, this data set is protected by the Creative Commons
#' License \url{https://creativecommons.org/licenses/by-nc-nd/2.5/}.
#' 
#' When publishing results obtained using this data set, the original authors
#' (Resnick et al, 1997) should be cited. In addition this package should be
#' cited as:
#' 
#' Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau, and
#' Martina Morris. 2003 \emph{statnet: Software tools for the Statistical
#' Modeling of Network Data} \cr \url{https://statnet.org}.
#' @seealso \code{\link[network]{network}},
#' \code{\link[network]{plot.network}}, \code{\link{ergm}},
#' \code{\link{faux.magnolia.high}}
#' @references
#' 
#' Hunter D.R., Goodreau S.M. and Handcock M.S. (2008). \emph{Goodness of Fit
#' of Social Network Models}, \emph{Journal of the American Statistical
#' Association}.
#' 
#' Resnick M.D., Bearman, P.S., Blum R.W. et al. (1997). \emph{Protecting
#' adolescents from harm. Findings from the National Longitudinal Study on
#' Adolescent Health}, \emph{Journal of the American Medical Association}, 278:
#' 823-32.
#' @source The data set is based upon a model fit to data from one school
#' community from the AddHealth Study, Wave I (Resnick et al., 1997). It was
#' constructed as follows:
#' 
#' A vector representing the sex of each student in the school was randomly
#' re-ordered.  The same was done with the students' response to questions on
#' race and grade.  These three attribute vectors were permuted independently.
#' Missing values for each were randomly assigned with weights determined by
#' the size of the attribute classes in the school.
#' 
#' The following \code{\link{ergm}} formula was used to fit a model to the
#' original data:
#' 
#' \preformatted{ ~ edges + nodefactor("Grade") + nodefactor("Race") +
#' nodefactor("Sex") + nodematch("Grade",diff=TRUE) +
#' nodematch("Race",diff=TRUE) + nodematch("Sex",diff=FALSE) +
#' gwdegree(1.0,fixed=TRUE) + gwesp(1.0,fixed=TRUE) + gwdsp(1.0,fixed=TRUE) }
#' 
#' The resulting model fit was then applied to a network with actors possessing
#' the permuted attributes and with the same number of edges as in the original
#' data.
#' 
#' The processes for handling missing data and defining the race attribute are
#' described in Hunter, Goodreau \& Handcock (2008).
#' @keywords datasets
NULL

#' Florentine Family Marriage and Business Ties Data as a "network" object
#' 
#' This is a data set of marriage and business ties among Renaissance
#' Florentine families. The data is originally from Padgett (1994) via
#' \code{UCINET} and stored as a \code{\link[network]{network}} object.
#' 
#' Breiger \& Pattison (1986), in their discussion of local role analysis, use
#' a subset of data on the social relations among Renaissance Florentine
#' families (person aggregates) collected by John Padgett from historical
#' documents. The two relations are business ties (\code{flobusiness} -
#' specifically, recorded financial ties such as loans, credits and joint
#' partnerships) and marriage alliances (\code{flomarriage}).
#' 
#' As Breiger \& Pattison point out, the original data are symmetrically coded.
#' This is acceptable perhaps for marital ties, but is unfortunate for the
#' financial ties (which are almost certainly directed). To remedy this, the
#' financial ties can be recoded as directed relations using some external
#' measure of power - for instance, a measure of wealth. Both graphs provide
#' vertex information on (1) \code{wealth} each family's net wealth in 1427 (in
#' thousands of lira); (2) \code{priorates} the number of priorates (seats on
#' the civic council) held between 1282- 1344; and (3) \code{totalties} the
#' total number of business or marriage ties in the total dataset of 116
#' families (see Breiger \& Pattison (1986), p 239).
#' 
#' Substantively, the data include families who were locked in a struggle for
#' political control of the city of Florence around 1430. Two factions were
#' dominant in this struggle: one revolved around the infamous Medicis (9), the
#' other around the powerful Strozzis (15).
#' 
#' @usage
#' data(florentine)
#' @docType data
#' @name florentine
#' @aliases flobusiness flomarriage
#' 
#' 
#' @seealso flo, network, plot.network, ergm
#' @references Wasserman, S. and Faust, K. (1994) \emph{Social Network
#' Analysis: Methods and Applications}, Cambridge University Press, Cambridge,
#' England.
#' 
#' Breiger R. and Pattison P. (1986). \emph{Cumulated social roles: The duality
#' of persons and their algebras}, Social Networks, 8, 215-256.
#' @source Padgett, John F. 1994. Marriage and Elite Structure in Renaissance
#' Florence, 1282-1500. Paper delivered to the Social Science History
#' Association.
#' @keywords datasets
NULL

#' Goodreau's four node network as a "network" object
#' 
#' This is an example thought of by Steve Goodreau. It is a directed network of
#' four nodes and five ties stored as a \code{\link[network]{network}} object.
#' 
#' It is interesting because the maximum likelihood estimator of the model with
#' out degree 3 in it exists, but the maximum psuedolikelihood estimator does
#' not.
#' 
#' @usage
#' data(g4)
#' @docType data
#' @name g4
#' 
#' @seealso florentine, network, plot.network, ergm
#' @source Steve Goodreau
#' @keywords datasets
#' @examples
#' 
#' data(g4)
#' summary(ergm(g4 ~ odegree(3), estimate="MPLE"))
#' summary(ergm(g4 ~ odegree(3), control=control.ergm(init=0)))
#' 
NULL


#' Kapferer's tailor shop data
#' 
#' This well-known social network dataset, collected by Bruce Kapferer in
#' Zambia from June 1965 to August 1965, involves interactions among workers in
#' a tailor shop as observed by Kapferer himself.
#'
#' An interaction is
#' defined by Kapferer as "continuous uninterrupted social activity involving
#' the participation of at least two persons"; only transactions that were
#' relatively frequent are recorded. All of the interactions in this particular
#' dataset are "sociational", as opposed to "instrumental".  Kapferer explains
#' the difference (p. 164) as follows:
#' 
#' "I have classed as transactions which were sociational in content those
#' where the activity was markedly convivial such as general conversation, the
#' sharing of gossip and the enjoyment of a drink together.  Examples of
#' instrumental transactions are the lending or giving of money, assistance at
#' times of personal crisis and help at work."
#' 
#' Kapferer also observed and recorded instrumental transactions, many of which
#' are unilateral (directed) rather than reciprocal (undirected), though those
#' transactions are not recorded here.  In addition, there was a second period
#' of data collection, from September 1965 to January 1966, but these data are
#' also not recorded here.  All data are given in Kapferer's 1972 book on pp.
#' 176-179.
#' 
#' During the first time period, there were 43 individuals working in this
#' particular tailor shop; however, the better-known dataset includes only
#' those 39 individuals who were present during both time collection periods.
#' (Missing are the workers named Lenard, Peter, Lazarus, and Laurent.) Thus,
#' we give two separate network datasets here: \code{kapferer} is the
#' well-known 39-individual dataset, whereas \code{kapferer2} is the full
#' 43-individual dataset.
#' 
#' 
#' @usage
#' data(kapferer)
#' @docType data
#' @name kapferer
#' @aliases kapferer kapferer2 tailor
#' @format Two \code{network} objects, \code{kapferer} and \code{kapferer2}.
#' The \code{kapferer} dataset contains only the 39 individuals who were
#' present at both data-collection time periods.  However, these data only
#' reflect data collected during the first period.  The individuals' names are
#' included as a nodal covariate called \code{names}.
#' @source Original source: Kapferer, Bruce (1972), Strategy and Transaction in
#' an African Factory, Manchester University Press.
#' @keywords datasets
NULL

#' Synthetic network with 20 nodes and 28 edges
#' 
#' This is a synthetic network of 20 nodes that is used as an example within
#' the \code{\link{ergm}} documentation. It has an interesting elongated shape
#' - reminencent of a chemical molecule.  It is stored as a
#' \code{\link[network]{network}} object.
#' 
#' 
#' @usage
#' data(molecule)
#' @docType data
#' @name molecule
#' @seealso florentine, sampson, network, plot.network, ergm
#' @keywords datasets
NULL

#' Longitudinal networks of positive affection within a monastery as a
#' "network" object
#' 
#' Three \code{\link{network}} objects containing the "liking" nominations of
#' Sampson's (1969) monks at the three time points.
#' 
#' Sampson (1969) recorded the social interactions among a group of monks while
#' he was a resident as an experimenter at the cloister.  During his stay, a
#' political "crisis in the cloister" resulted in the expulsion of four
#' monks-- namely, the three "outcasts," Brothers Elias, Simplicius, Basil, and
#' the leader of the "young Turks," Brother Gregory.  Not long after Brother
#' Gregory departed, all but one of the "young Turks" left voluntarily:
#' Brothers John Bosco, Albert, Boniface, Hugh, and Mark.  Then, all three of
#' the "waverers" also left: First, Brothers Amand and Victor, then later
#' Brother Romuald.  Eventually, Brother Peter and Brother Winfrid also left,
#' leaving only four of the original group.
#' 
#' Of particular interest are the data on positive affect relations
#' ("liking," using the terminology later adopted by White et al. (1976)), in
#' which each monk was asked if he had positive relations to each of the other
#' monks. Each monk ranked only his top three choices (or four, in the case of
#' ties) on "liking".  Here, we consider a directed edge from monk A to monk
#' B to exist if A nominated B among these top choices.
#' 
#' The data were gathered at three times to capture changes in group sentiment
#' over time. They represent three time points in the period during which a new
#' cohort had entered the monastery near the end of the study but before the
#' major conflict began.  These three time points are labeled T2, T3, and T4 in
#' Tables D5 through D16 in the appendices of Sampson's 1969 dissertation.  and
#' the corresponding network data sets are named \code{samplk1},
#' \code{samplk2}, and \code{samplk3}, respectively.
#' 
#' See also the data set \code{\link{sampson}} containing the time-aggregated
#' graph \code{samplike}.
#' 
#' \code{samplk3} is a data set of Hoff, Raftery and Handcock (2002).
#' 
#' The data sets are stored as \code{\link[network]{network}} objects with
#' three vertex attributes:
#' 
#' \describe{ \item{group}{Groups of novices as classified by Sampson, that is,
#' "Loyal", "Outcasts", and "Turks", but with a fourth group called the
#' "Waverers" by White et al. (1975) that comprises two of the original Loyal
#' opposition and one of the original Outcasts. See the \code{\link{samplike}}
#' data set for the original classifications of these three waverers.}
#' \item{cloisterville}{An indicator of attendance in the minor seminary of
#' "Cloisterville" before coming to the monastery.} \item{vertex.names}{The
#' given names of the novices. NB: These names have been corrected as of
#' \code{ergm} version 3.6.1.} } This data set is standard in the social
#' network analysis literature, having been modeled by Holland and Leinhardt
#' (1981), Reitz (1982), Holland, Laskey and Leinhardt (1983), Fienberg, Meyer,
#' and Wasserman (1981), and Hoff, Raftery, and Handcock (2002), among others.
#' This is only a small piece of the data collected by Sampson.
#' 
#' This data set was updated for version 2.5 (March 2012) to add the
#' \code{cloisterville} variable and refine the names. This information is from
#' de Nooy, Mrvar, and Batagelj (2005). The original vertex names were:
#' Romul_10, Bonaven_5, Ambrose_9, Berth_6, Peter_4, Louis_11, Victor_8,
#' Winf_12, John_1, Greg_2, Hugh_14, Boni_15, Mark_7, Albert_16, Amand_13,
#' Basil_3, Elias_17, Simp_18. The numbers indicate the ordering used in the
#' original dissertation of Sampson (1969).
#' 
#' @usage
#' data(samplk)
#' @docType data
#' @name samplk
#' @aliases samplk samplk1 samplk2 samplk3
#' @section Mislabeling in Versions Prior to 3.6.1: In \code{ergm} versions
#' 3.6.0 and earlier, The adjacency matrices of the \code{\link{samplike}},
#' \code{\link{samplk1}}, \code{\link{samplk2}}, and \code{\link{samplk3}}
#' networks reflected the original Sampson (1969) ordering of the names even
#' though the vertex labels used the name order of de Nooy, Mrvar, and Batagelj
#' (2005). That is, in \code{ergm} version 3.6.0 and earlier, the vertices were
#' mislabeled. The correct order is the same one given in Tables D5, D9, and
#' D13 of Sampson (1969): John Bosco, Gregory, Basil, Peter, Bonaventure,
#' Berthold, Mark, Victor, Ambrose, Romauld (Sampson uses both spellings
#' "Romauld" and "Ramauld" in the dissertation), Louis, Winfrid, Amand, Hugh,
#' Boniface, Albert, Elias, Simplicius. By contrast, the order given in
#' \code{ergm} version 3.6.0 and earlier is: Ramuald, Bonaventure, Ambrose,
#' Berthold, Peter, Louis, Victor, Winfrid, John Bosco, Gregory, Hugh,
#' Boniface, Mark, Albert, Amand, Basil, Elias, Simplicius.
#' @seealso sampson, florentine, network, plot.network, ergm
#' @references White, H.C., Boorman, S.A. and Breiger, R.L. (1976).
#' \emph{Social structure from multiple networks. I. Blockmodels of roles and
#' positions.} American Journal of Sociology, 81(4), 730-780.
#' 
#' Wouter de Nooy, Andrej Mrvar, Vladimir Batagelj (2005) \emph{Exploratory
#' Social Network Analysis with Pajek}, Cambridge: Cambridge University Press
#' @source Sampson, S.~F. (1968), \emph{A novitiate in a period of change: An
#' experimental and case study of relationships,} Unpublished Ph.D.
#' dissertation, Department of Sociology, Cornell University.
#' 
#' \url{http://vlado.fmf.uni-lj.si/pub/networks/data/esna/sampson.htm}
#' @keywords datasets
NULL

#' Cumulative network of positive affection within a monastery as a "network"
#' object
#' 
#' A \code{\link{network}} object containing the cumulative "liking"
#' nominations of Sampson's (1969) monks over the three time points.
#' 
#' Sampson (1969) recorded the social interactions among a group of monks while
#' he was a resident as an experimenter at the cloister.  During his stay, a
#' political "crisis in the cloister" resulted in the expulsion of four
#' monks-- namely, the three "outcasts," Brothers Elias, Simplicius, Basil, and
#' the leader of the "young Turks," Brother Gregory.  Not long after Brother
#' Gregory departed, all but one of the "young Turks" left voluntarily:
#' Brothers John Bosco, Albert, Boniface, Hugh, and Mark.  Then, all three of
#' the "waverers" also left: First, Brothers Amand and Victor, then later
#' Brother Romuald.  Eventually, Brother Peter and Brother Winfrid also left,
#' leaving only four of the original group.
#' 
#' Of particular interest are the data on positive affect relations
#' ("liking," using the terminology later adopted by White et al. (1976)), in
#' which each monk was asked if he had positive relations to each of the other
#' monks. Each monk ranked only his top three choices (or four, in the case of
#' ties) on "liking".  Here, we consider a directed edge from monk A to monk
#' B to exist if A nominated B among these top choices.
#' 
#' The data were gathered at three times to capture changes in group sentiment
#' over time. They represent three time points in the period during which a new
#' cohort had entered the monastery near the end of the study but before the
#' major conflict began.  These three time points are labeled T2, T3, and T4 in
#' Tables D5 through D16 in the appendices of Sampson's 1969 dissertation.  The
#' \code{samplike} data set is the time-aggregated network.  Thus, a tie from
#' monk A to monk B exists if A nominated B as one of his three (or four, in
#' case of ties) best friends at any of the three time points.
#' 
#' See also the data sets \code{\link{samplk1}}, \code{\link{samplk2}}, and
#' \code{\link{samplk3}}, containing the networks at each of the three
#' individual time points.
#' 
#' The data set is stored as a \code{\link[network]{network}} object with three
#' vertex attributes:
#' 
#' \describe{ \item{group}{Groups of novices as classified by Sampson:
#' "Loyal", "Outcasts", and "Turks".} \item{cloisterville}{An indicator
#' of attendance in the minor seminary of "Cloisterville" before coming to
#' the monastery.} \item{vertex.names}{The given names of the novices.  NB:
#' These names have been corrected as of \code{ergm} version 3.6.1; see details
#' below.} } In addition, the data set has an edge attribute,
#' \code{nominations}, giving the number of times (out of 3) that monk A
#' nominated monk B.
#' 
#' This data set is standard in the social network analysis literature, having
#' been modeled by Holland and Leinhardt (1981), Reitz (1982), Holland, Laskey
#' and Leinhardt (1983), Fienberg, Meyer, and Wasserman (1981), and Hoff,
#' Raftery, and Handcock (2002), among others. This is only a small piece of
#' the data collected by Sampson.
#' 
#' This data set was updated for version 2.5 (March 2012) to add the
#' \code{cloisterville} variable and refine the names. This information is from
#' de Nooy, Mrvar, and Batagelj (2005). The original vertex names were:
#' Romul_10, Bonaven_5, Ambrose_9, Berth_6, Peter_4, Louis_11, Victor_8,
#' Winf_12, John_1, Greg_2, Hugh_14, Boni_15, Mark_7, Albert_16, Amand_13,
#' Basil_3, Elias_17, Simp_18. The numbers indicate the ordering used in the
#' original dissertation of Sampson (1969).
#' 
#' @usage
#' data(sampson)
#' @docType data
#' @name sampson
#' @aliases sampson samplike
#' @section Mislabeling in Versions Prior to 3.6.1: In \code{ergm} version
#' 3.6.0 and earlier, The adjacency matrices of the \code{\link{samplike}},
#' \code{\link{samplk1}}, \code{\link{samplk2}}, and \code{\link{samplk3}}
#' networks reflected the original Sampson (1969) ordering of the names even
#' though the vertex labels used the name order of de Nooy, Mrvar, and Batagelj
#' (2005). That is, in \code{ergm} version 3.6.0 and earlier, the vertices were
#' mislabeled. The correct order is the same one given in Tables D5, D9, and
#' D13 of Sampson (1969): John Bosco, Gregory, Basil, Peter, Bonaventure,
#' Berthold, Mark, Victor, Ambrose, Romauld (Sampson uses both spellings
#' "Romauld" and "Ramauld" in the dissertation), Louis, Winfrid, Amand, Hugh,
#' Boniface, Albert, Elias, Simplicius. By contrast, the order given in
#' \code{ergm} version 3.6.0 and earlier is: Ramuald, Bonaventure, Ambrose,
#' Berthold, Peter, Louis, Victor, Winfrid, John Bosco, Gregory, Hugh,
#' Boniface, Mark, Albert, Amand, Basil, Elias, Simplicius.
#' @seealso florentine, network, plot.network, ergm
#' @references White, H.C., Boorman, S.A. and Breiger, R.L. (1976).
#' \emph{Social structure from multiple networks. I. Blockmodels of roles and
#' positions.} American Journal of Sociology, 81(4), 730-780.
#' 
#' Wouter de Nooy, Andrej Mrvar, Vladimir Batagelj (2005) \emph{Exploratory
#' Social Network Analysis with Pajek}, Cambridge: Cambridge University Press
#' @source Sampson, S.~F. (1968), \emph{A novitiate in a period of change: An
#' experimental and case study of relationships,} Unpublished Ph.D.
#' dissertation, Department of Sociology, Cornell University.
#' 
#' \url{http://vlado.fmf.uni-lj.si/pub/networks/data/esna/sampson.htm}
#' @keywords datasets
NULL
