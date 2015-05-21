#ifndef WTCD_H
#define WTCD_H

#include "wtMCMC.h"

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */
void WtCD_wrapper(int *dnumnets, int *nedges, 
		    int *tails, int *heads, double *weights,
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHproposaltype, char **MHproposalpackage,
		    double *inputs, double *theta0, int *samplesize, int *CDparams,
		    double *sample,
		    int *fVerbose, 
		    int *status);
WtMCMCStatus WtCDSample(WtMHproposal *MHp,
			double *theta, double *networkstatistics, 
			int samplesize, int *CDparams, Vertex *undotail, Vertex *undohead, double *undoweight,
			int fVerbose,
			WtNetwork *nwp, WtModel *m, double *extraworkspace);
WtMCMCStatus WtCDStep(WtMHproposal *MHp,
		      double *theta, double *statistics, 
		      int *CDparams, int *staken, Vertex *undotail, Vertex *undohead, double *undoweight,
		      int fVerbose,
		      WtNetwork *nwp, WtModel *m, double *extraworkspace);

#endif
