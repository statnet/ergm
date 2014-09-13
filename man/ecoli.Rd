%  File man/ecoli.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ecoli}
\docType{data}
\title{Two versions of an E. Coli network dataset}
\alias{ecoli}
\alias{ecoli1}
\alias{ecoli2}
\description{
This network data set comprises two versions of a biological
network in which the nodes are operons in \emph{Escherichia Coli}
and a directed edge from one node to another indicates that 
the first encodes the transcription factor that regulates the second.
}

\usage{data(ecoli)}

\details{
The network object \code{ecoli1} is directed, with 423 nodes and 519 arcs.
The object \code{ecoli2} is an undirected version of the same network, in 
which all arcs are treated as edges and the five isolated nodes (which exhibit
only self-regulation in \code{ecoli1}) are removed, leaving 418 nodes.
}

\source{
The data set is based on the RegulonDB network (Salgado et al, 2001)
and was modified by Shen-Orr et al (2002).  

}


\references{

Salgado et al (2001), Regulondb (version 3.2):  Transcriptional Regulation
and Operon Organization in Escherichia Coli K-12,
\emph{Nucleic Acids Research}, 29(1): 72-74.

Shen-Orr et al (2002), Network Motifs in the Transcriptional Regulation
Network of Escerichia Coli, \emph{Nature Genetics},
31(1): 64-68.

%Saul and Filkov (2007)

%Hummel et al (2010)

}

\keyword{datasets}
\section{Licenses and Citation}{
When publishing results obtained using this data set, the original authors 
(Salgado et al, 2001; Shen-Orr et al, 2002) should
be cited, along with this \code{R} package.
}
