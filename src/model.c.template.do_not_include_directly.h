/*  File src/model.c.template.do_not_include_directly.h in package ergm, part
 *  of the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include <string.h>
#include "ergm_omp.h"
#include "ergm_util.h"

void ETYPE(On, NetworkEdgeChangeUWrap)(Vertex tail, Vertex head, IFEWT(EWTTYPE weight,) void *mtp, ETYPE(Network) *nwp, EWTTYPE edgestate){
  ((ETYPE(ModelTerm) *) mtp)->u_func(tail, head, IFEWT(weight,) mtp, nwp, edgestate);
}

/*
  ETYPE(InitStats)
  A helper's helper function to initialize storage for functions that use it.
*/
static inline void ETYPE(InitStats)(ETYPE(Network) *nwp, ETYPE(Model) *m){

  /* This function must do things in very specific order:

   1. Since dependent terms go before the dependency, the
      initialization must be performed in reverse order.

   2. Since a dependt term relies on the pre-toggle state of the
      dependency, the updating must be performed in forward order.

   3. Since an i_function can choose to delete its own u_function, we
      can't add callbacks until after the i_functions have been
      called.

   4. Some terms may have subterms that add callbacks to the same
      network, and the subterms' u_functions must be called *after*
      the u_functions of the terms that depend on them (per rule 2).

   Therefore, this code first stores the current position in the
   callback list, then initializes the terms, then adds callbacks in
   front of the new callbacks, i.e., at the then-current
   position. Repeatedly adding the terms at that position in reverse
   order is slightly inefficient, but it has to be done infrequently.

   */
  unsigned int on_edge_change_pos = nwp->n_on_edge_change; // Save the current position.

  // Iterate in reverse, so that auxliary terms get initialized first.
  ETYPE(EXEC_THROUGH_TERMS_INREVERSE)(m, {
      if(!m->noinit_s || !mtp->s_func){ // Skip if noinit_s is set and s_func is present.
        double *dstats = mtp->dstats;
        mtp->dstats = NULL; // Trigger segfault if i_func tries to write to change statistics.
        if(mtp->i_func)
          (*(mtp->i_func))(mtp, nwp);  /* Call i_??? function */
        else if(mtp->u_func) /* No initializer but an updater -> uses a 1-function implementation. */
          (*(mtp->u_func))(0, 0, IFEWT(0,) mtp, nwp, 0);  /* Call u_??? function */
        mtp->dstats = dstats;
      }
      // Now, bind the term to the network through the callback API.
      if(mtp->u_func && (!m->noinit_s || !mtp->s_func)) // Skip if noinit_s is set and s_func is present.
        ETYPE(AddOn, NetworkEdgeChange)(nwp, ETYPE(On, NetworkEdgeChangeUWrap), mtp, on_edge_change_pos);
    });
}

/*
  ETYPE(DestroyStats)
  A helper's helper function to finalize storage for functions that use it.
*/
static inline void ETYPE(DestroyStats)(ETYPE(Network) *nwp, ETYPE(Model) *m){
  unsigned int i=0;
  ETYPE(EXEC_THROUGH_TERMS)(m, {
      if(!m->noinit_s || !mtp->s_func){ // Skip if noinit_s is set and s_func is present.
        if(mtp->u_func)
          ETYPE(DeleteOn, NetworkEdgeChange)(nwp, ETYPE(On, NetworkEdgeChangeUWrap), mtp);
        if(mtp->f_func)
          (*(mtp->f_func))(mtp, nwp);  /* Call f_??? function */
      }
      R_Free(m->dstatarray[i]);
      R_Free(mtp->statcache);
      if(mtp->storage){
        R_Free(mtp->storage);
        mtp->storage = NULL;
      }
      i++;
    });
}

/*****************
  void ETYPE(ModelDestroy)
******************/
void ETYPE(ModelDestroy)(ETYPE(Network) *nwp, ETYPE(Model) *m)
{
  ETYPE(DestroyStats)(nwp, m);

  for(unsigned int i=0; i < m->n_aux; i++)
    if(m->termarray[0].aux_storage[i]!=NULL){
      R_Free(m->termarray[0].aux_storage[i]);
      m->termarray[0].aux_storage[i] = NULL;
  }

  if(m->n_terms && m->termarray[0].aux_storage!=NULL){
    R_Free(m->termarray[0].aux_storage);
  }

  ETYPE(EXEC_THROUGH_TERMS)(m, {
      if(mtp->aux_storage!=NULL)
	mtp->aux_storage=NULL;
    });

  R_Free(m->dstatarray);
  R_Free(m->termarray);
  R_Free(m->workspace_backup);
  R_Free(m);
}

/*****************
 int ETYPE(ModelInitialize)

 Allocate and initialize the ETYPE(ModelTerm) structures, each of which contains
 all necessary information about how to compute one term in the model.
*****************/
ETYPE(Model)* ETYPE(ModelInitialize) (SEXP mR, SEXP ext_state, ETYPE(Network) *nwp, Rboolean noinit_s){
  SEXP terms = getListElement(mR, "terms");
  if(ext_state == R_NilValue) ext_state = NULL;

  ETYPE(Model) *m = (ETYPE(Model) *) R_Calloc(1, ETYPE(Model));
  unsigned int n_terms = m->n_terms = length(terms);
  m->termarray = (ETYPE(ModelTerm) *) R_Calloc(n_terms, ETYPE(ModelTerm));
  m->dstatarray = (double **) R_Calloc(n_terms, double *);
  m->n_stats = 0;
  m->n_aux = 0;
  m->n_u = 0;
  m->noinit_s = noinit_s;
  m->R = mR;
  m->ext_state = ext_state;
  for (unsigned int l=0; l < m->n_terms; l++) {
    ETYPE(ModelTerm) *thisterm = m->termarray + l;
    thisterm->R = VECTOR_ELT(terms, l);

      /* Initialize storage and term functions to NULL. */
      thisterm->storage = NULL;
      thisterm->aux_storage = NULL;
      thisterm->ext_state = NULL;
      thisterm->d_func = NULL;
      thisterm->c_func = NULL;
      thisterm->s_func = NULL;
      thisterm->i_func = NULL;
      thisterm->u_func = NULL;
      thisterm->f_func = NULL;
      thisterm->w_func = NULL;
      thisterm->x_func = NULL;

      /* First, obtain the term name and library: fnames points to a
      single character string, consisting of the names of the selected
      options concatenated together and separated by spaces.  This is
      passed by the calling R function.  These names are matched with
      their respective C functions that calculate the appropriate
      statistics.  Similarly, sonames points to a character string
      containing the names of the shared object files associated with
      the respective functions.*/
      const char *fname = FIRSTCHAR(getListElement(thisterm->R, "name")),
        *sn = FIRSTCHAR(getListElement(thisterm->R, "pkgname"));

      /* Extract the required string information from the relevant sources */
      char *fn = R_Calloc(strlen(fname)+3, char);
      fn[1]='_';
      strcpy(fn+2, fname);
      /* fn is now the string ' _[name]', where [name] is fname */

      /* Extract the term inputs. */

      /* Double input vector with an optional attribute shift. */
      SEXP tmp = getListElement(thisterm->R, "inputs");
      thisterm->ninputparams = length(tmp);
      thisterm->inputparams = thisterm->ninputparams ? REAL(tmp) : NULL;

      tmp = getAttrib(tmp, install("ParamsBeforeCov"));
      unsigned int offset = length(tmp) ? asInteger(tmp): 0;  /* Set offset for attr vector */
      thisterm->attrib = thisterm->ninputparams ? thisterm->inputparams + offset : NULL; /* Ptr to attributes */

      /* Integer input vector with an optional attribute shift. */
      tmp = getListElement(thisterm->R, "iinputs");
      thisterm->niinputparams = length(tmp);
      thisterm->iinputparams = thisterm->niinputparams ? INTEGER(tmp) : NULL;

      tmp = getAttrib(tmp, install("ParamsBeforeCov"));
      offset = length(tmp) ? asInteger(tmp): 0;  /* Set offset for attr vector */
      thisterm->iattrib = thisterm->niinputparams ? thisterm->iinputparams + offset : NULL; /* Ptr to attributes */

      /* Number of statistics. */
      thisterm->nstats = length(getListElement(thisterm->R, "coef.names")); /* If >0, # of statistics returned. If ==0 an auxiliary statistic. */

      /* Set auxiliary counts and values. */
      tmp = getAttrib(thisterm->R, install("aux.slots"));
      thisterm->n_aux = length(tmp);
      thisterm->aux_slots = (unsigned int *) INTEGER(tmp);

      /* Empty network statistics. */
      tmp = getListElement(thisterm->R, "emptynwstats");
      thisterm->emptynwstats = isNULL(tmp) ? NULL : REAL(tmp);

      /*  Update the running total of statistics */
      m->n_stats += thisterm->nstats;
      m->dstatarray[l] = (double *) R_Calloc(thisterm->nstats, double);
      thisterm->dstats = m->dstatarray[l];  /* This line is important for
                                               eventually freeing up allocated
					       memory, since thisterm->dstats
					       can be modified but
					       m->dstatarray[l] cannot be.  */
      thisterm->statcache = (double *) R_Calloc(thisterm->nstats, double);

      if(ext_state) thisterm->ext_state = VECTOR_ELT(ext_state, l);

      /* If the term's nstats==0, it is auxiliary: it does not affect
	 acceptance probabilities or contribute any
	 statistics. Therefore, its d_ and s_ functions are never
	 called and are not initialized. It is only used for its u_
	 function. Therefore, the following code is only run when
	 thisterm->nstats>0. */
      if(thisterm->nstats){
	/*  Most important part of the ETYPE(ModelTerm):  A pointer to a
	    function that will compute the change in the network statistic of
	    interest for a particular edge toggle.  This function is obtained by
	    searching for symbols associated with the object file with prefix
	    sn, having the name fn.  Assuming that one is found, we're golden.*/
	fn[0]='c';
	thisterm->c_func =
	  (void (*)(Vertex, Vertex, IFEWT(EWTTYPE,) ETYPE(ModelTerm)*, ETYPE(Network)*, EWTTYPE))
	  R_FindSymbol(fn,sn,NULL);

        fn[0]='d';
        thisterm->d_func =
          (void (*)(Edge, Vertex*, Vertex*, IFEWT(EWTTYPE*,) ETYPE(ModelTerm)*, ETYPE(Network)*))
          R_FindSymbol(fn,sn,NULL);

	/* Optional function to compute the statistic of interest for
	   the network given. It can be more efficient than going one
	   edge at a time. */
	fn[0]='s';
	thisterm->s_func =
	  (void (*)(ETYPE(ModelTerm)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);

        if(thisterm->c_func==NULL && thisterm->d_func==NULL && thisterm->s_func==NULL){
          error("Error in C model initialization: term with functions %s::%s is declared to have statistics but does not appear to have a change, a difference, or a summary function. Memory has not been deallocated, so restart R sometime soon.\n",sn,fn+2);
	}

	/* Optional function to compute the statistic of interest for
	   the empty network (over and above the constant value if
	   given) and taking into account the extended state. */
	fn[0]='z';
	thisterm->z_func =
	  (void (*)(ETYPE(ModelTerm)*, ETYPE(Network)*, Rboolean)) R_FindSymbol(fn,sn,NULL);
      }else m->n_aux++;

      /* Optional functions to store persistent information about the
	 network state between calls to d_ functions. */
      fn[0]='u';
      if((thisterm->u_func =
	  (void (*)(Vertex, Vertex, IFEWT(EWTTYPE,) ETYPE(ModelTerm)*, ETYPE(Network)*, EWTTYPE)) R_FindSymbol(fn,sn,NULL))!=NULL) m->n_u++;

      /* Optional-optional functions to initialize and finalize the
	 term's storage, and the "eXtension" function to allow an
	 arbitrary "signal" to be sent to a statistic. */

      fn[0]='i';
      thisterm->i_func =
	(void (*)(ETYPE(ModelTerm)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);

      fn[0]='f';
      thisterm->f_func =
	(void (*)(ETYPE(ModelTerm)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);

      /* If it's an auxiliary, then it needs an i_function or a
	 u_function, or it's not doing anything. */
      if(thisterm->nstats==0 && (thisterm->i_func==NULL && thisterm->u_func==NULL)){
          error("Error in C model initialization: term with functions %s::%s is declared to have no statistics but does not appear to have an updater function, so does not do anything. Memory has not been deallocated, so restart R sometime soon.\n",sn,fn+2);
      }

      fn[0]='w';
      thisterm->w_func =
	(SEXP (*)(ETYPE(ModelTerm)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);

      fn[0]='x';
      thisterm->x_func =
	(void (*)(unsigned int type, void *data, ETYPE(ModelTerm)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);

      if(!ext_state && (thisterm->w_func)) error("Error in ModelInitialize: not provided with extended state, but model terms with functions %s::%s requires extended state. This should normally be caught sooner. This limitation may be removed in the future.  Memory has not been deallocated, so restart R sometime soon.\n",sn,fn+2);

      /*Clean up by freeing fn*/
      R_Free(fn);
  }

  m->workspace_backup = m->workspace = (double *) R_Calloc(m->n_stats, double);

  unsigned int pos = 0;
  ETYPE(FOR_EACH_TERM)(m){
    mtp->statspos = pos;
    pos += mtp->nstats;
  }

  /* Allocate auxiliary storage and put a pointer to it on every model term. */
  if(m->n_aux){
    m->termarray[0].aux_storage = (void *) R_Calloc(m->n_aux, void *);
    for(unsigned int l=1; l < n_terms; l++)
      m->termarray[l].aux_storage = m->termarray[0].aux_storage;
  }

  /* Trigger initial storage update */
  ETYPE(InitStats)(nwp, m);

  /* Now, check that no term exports both a d_ and a c_
     function. TODO: provide an informative "traceback" to which term
     caused the problem.*/
  ETYPE(FOR_EACH_TERM)(m){
    if(mtp->c_func && mtp->d_func) error("A term exports both a change and a difference function.  Memory has not been deallocated, so restart R sometime soon.\n");
  }

  return m;
}

void ETYPE(ChangeStatsDo)(unsigned int ntoggles, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,)
                   ETYPE(Network) *nwp, ETYPE(Model) *m){
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */

  /* Make a pass through terms with d_functions. */
  ETYPE(EXEC_THROUGH_TERMS_INTO)(m, m->workspace, {
      mtp->dstats = dstats; /* Stuck the change statistic here.*/
      if(mtp->c_func==NULL && mtp->d_func)
	(*(mtp->d_func))(ntoggles, tails, heads, IFEWT(weights,)
			 mtp, nwp);  /* Call d_??? function */
    });
  /* Notice that mtp->dstats now points to the appropriate location in
     m->workspace. */

  /* Put the original destination dstats back unless there's only one
     toggle. */
  if(ntoggles!=1){
    unsigned int i = 0;
    ETYPE(EXEC_THROUGH_TERMS)(m, {
	mtp->dstats = m->dstatarray[i];
	i++;
      });
  }

  /* Make a pass through terms with c_functions. */
  IFELSEEWT(FOR_EACH_TOGGLE,
            int toggle;
            FOR_EACH_TOGGLE(toggle)){
    IFELSEEWT(GETTOGGLEINFO(), Rboolean edgestate = IS_OUTEDGE(tails[toggle], heads[toggle]));

    ergm_PARALLEL_FOR_LIMIT(m->n_terms)
    ETYPE(EXEC_THROUGH_TERMS_INTO)(m, m->workspace, {
	if(mtp->c_func){
	  if(ntoggles!=1) ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))IFELSEEWT((TAIL, HEAD, NEWWT,
                                     mtp, nwp, OLDWT),
                                    (tails[toggle], heads[toggle],
                                     mtp, nwp, edgestate));  /* Call d_??? function */
	  if(ntoggles!=1){
            addonto(dstats, mtp->dstats, N_CHANGE_STATS);
	  }
	}
      });

    /* Update storage and network */
    IFELSEEWT(IF_MORE_TO_COME{SETWT_WITH_BACKUP();},
             IF_MORE_TO_COME(toggle){TOGGLE_KNOWN(tails[toggle],heads[toggle], edgestate);})

  }
}


void ETYPE(ChangeStatsUndo)(unsigned int ntoggles, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,)
                       ETYPE(Network) *nwp, ETYPE(Model) *m){
  IFELSEEWT(
  UNDO_PREVIOUS{
    GETOLDTOGGLEINFO();
    SETWT(TAIL,HEAD,weights[TOGGLEIND]);
    weights[TOGGLEIND]=OLDWT;
  },
    int toggle = ntoggles;
  UNDO_PREVIOUS(toggle){
    TOGGLE(tails[toggle],heads[toggle]);
  })

}


/*
  ETYPE(ChangeStats)
  A helper's helper function to compute change statistics.
  The vector of changes is written to m->workspace.
*/
void ETYPE(ChangeStats)(unsigned int ntoggles, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,)
				 ETYPE(Network) *nwp, ETYPE(Model) *m){
  ETYPE(ChangeStatsDo)(ntoggles, tails, heads, IFEWT(weights,) nwp, m);
  ETYPE(ChangeStatsUndo)(ntoggles, tails, heads, IFEWT(weights,) nwp, m);
}

/*
  ETYPE(ChangeStats1)
  A simplified version of ETYPE(ChangeStats) for exactly one change.
*/
void ETYPE(ChangeStats1)(Vertex tail, Vertex head, IFEWT(EWTTYPE weight,)
                    ETYPE(Network) *nwp, ETYPE(Model) *m, EWTTYPE edgestate){
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */

  /* Make a pass through terms with c_functions. */
  ergm_PARALLEL_FOR_LIMIT(m->n_terms)
    ETYPE(EXEC_THROUGH_TERMS_INTO)(m, m->workspace, {
        mtp->dstats = dstats; /* Stuck the change statistic here.*/
        if(mtp->c_func){
          (*(mtp->c_func))(tail, head, IFEWT(weight,)
                           mtp, nwp, edgestate);  /* Call c_??? function */
        }else if(mtp->d_func){
          (*(mtp->d_func))(1, &tail, &head, IFEWT(&weight,)
                           mtp, nwp);  /* Call d_??? function */
        }
      });
}


/*
  ETYPE(ZStats)
  Call baseline statistics calculation (for extended state).
*/
void ETYPE(ZStats)(ETYPE(Network) *nwp, ETYPE(Model) *m, Rboolean skip_s){
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */

  /* Make a pass through terms with c_functions. */
  ergm_PARALLEL_FOR_LIMIT(m->n_terms)
    ETYPE(EXEC_THROUGH_TERMS_INTO)(m, m->workspace, {
        mtp->dstats = dstats; /* Stuck the change statistic here.*/
        if(!skip_s || mtp->s_func==NULL)
          if(mtp->z_func)
            (*(mtp->z_func))(mtp, nwp, skip_s);  /* Call z_??? function */
      });
}

/*
  ETYPE(EmptyNetworkStats)
  Extract constant empty network stats.
*/
void ETYPE(EmptyNetworkStats)(ETYPE(Model) *m, Rboolean skip_s){
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */

  ETYPE(EXEC_THROUGH_TERMS_INTO)(m, m->workspace, {
      if(!skip_s || mtp->s_func==NULL){
        if(mtp->emptynwstats)
          memcpy(dstats, mtp->emptynwstats, mtp->nstats*sizeof(double));
      }});
}

/****************
 void ETYPE(SummStats)

Compute summary statistics. It has two modes:
* nwp is empty and m is initialized consistently with nwp -> use edgelist
* nwp is not empty n_edges=0, and m does not have to be initialized consistently with nwp -> use nwp (making temporary copies of nwp and reinitializing m)
*****************/
void ETYPE(SummStats)(Edge n_edges, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) ETYPE(Network) *nwp, ETYPE(Model) *m){
  Rboolean mynet;
  double *stats;
  if(EDGECOUNT(nwp)){
    if(n_edges) error("SummStats must be passed either an empty network and a list of edges or a non-empty network and no edges.");
    /* The following code is pretty inefficient, but it'll do for now. */
    /* Grab network state and output workspace. */
    n_edges = EDGECOUNT(nwp);
    /* Use R's memory management to make the routine interruptible.

       TODO: Check how much overhead this incurs over and above
       in-house on.exit() memory management.
    */
    tails = (Vertex *) INTEGER(PROTECT(allocVector(INTSXP, n_edges)));
    heads = (Vertex *) INTEGER(PROTECT(allocVector(INTSXP, n_edges)));
    IFEWT(weights = REAL(PROTECT(allocVector(REALSXP, n_edges))));

    ETYPE(EdgeTree2EdgeList)(tails, heads, IFEWT(weights,) nwp, n_edges);
    stats = m->workspace;

    /* Replace the model and network with an empty one. */
    nwp = ETYPE(NetworkInitializeLike)(nwp);
    m = ETYPE(ModelInitialize)(m->R, m->ext_state, nwp, TRUE);
    mynet = TRUE;
  }else{
    /* Use R's memory management to make the routine interruptible.

       TODO: Check how much overhead this incurs over and above
       in-house on.exit() memory management.
    */
    stats = REAL(PROTECT(allocVector(REALSXP, m->n_stats)));
    mynet = FALSE;
  }

  memset(stats, 0, m->n_stats*sizeof(double));

  ETYPE(EmptyNetworkStats)(m, TRUE);
  addonto(stats, m->workspace, m->n_stats);
  ETYPE(ZStats)(nwp, m, TRUE);
  addonto(stats, m->workspace, m->n_stats);

  ETYPE(DetShuffleEdges)(tails,heads,IFEWT(weights,)n_edges); /* Shuffle edgelist. */

  Edge ntoggles = n_edges; // So that we can use the macros

  /* Calculate statistics for terms that don't have c_functions or s_functions.  */
  ETYPE(EXEC_THROUGH_TERMS_INTO)(m, stats, {
      if(mtp->s_func==NULL && mtp->c_func==NULL && mtp->d_func){
	(*(mtp->d_func))(ntoggles, tails, heads, IFEWT(weights,)
			 mtp, nwp);  /* Call d_??? function */
        addonto(dstats, mtp->dstats, N_CHANGE_STATS);
      }
    });

  /* Calculate statistics for terms that have c_functions but not s_functions.  */
  IFELSEEWT(FOR_EACH_TOGGLE, for(Edge e=0; e<n_edges; e++)){
    IFELSEEWT(GETNEWTOGGLEINFO(), Vertex t=TAIL(e); Vertex h=HEAD(e));

    ergm_PARALLEL_FOR_LIMIT(m->n_terms)
    ETYPE(EXEC_THROUGH_TERMS_INTO)(m, stats, {
	if(mtp->s_func==NULL && mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))IFELSEEWT((TAIL, HEAD, NEWWT,
                                     mtp, nwp, 0),
                                    (t, h,
                                     mtp, nwp, FALSE));  /* Call c_??? function */
          addonto(dstats, mtp->dstats, N_CHANGE_STATS);
	}
      });

    /* Update storage and network */
    IFELSEEWT(SETWT(TAIL, HEAD, NEWWT),
              TOGGLE_KNOWN(t, h, FALSE));
  }

  /* Calculate statistics for terms have s_functions  */
  ETYPE(EXEC_THROUGH_TERMS_INTO)(m, stats, {
      if(mtp->s_func){
	ZERO_ALL_CHANGESTATS();
	(*(mtp->s_func))(mtp, nwp);  /* Call d_??? function */
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	  dstats[k] = mtp->dstats[k]; // Overwrite, not accumulate.
	}
      }
    });

  if(mynet){
    ETYPE(ModelDestroy)(nwp,m);
    ETYPE(NetworkDestroy)(nwp);
    UNPROTECT(IFELSEEWT(3,2));
  }else{
    ETYPE(DetUnShuffleEdges)(tails,heads,IFEWT(weights,)n_edges); /* Unshuffle edgelist. */
    memcpy(m->workspace, stats, m->n_stats*sizeof(double));
    UNPROTECT(1);
  }
}


/****************
 void ETYPE(SummStatsS)

Compute summary statistics, but only for terms that only have s_ statistics. Set the rest to NA_REAL.
*****************/
void ETYPE(SummStatsS)(ETYPE(Network) *nwp, ETYPE(Model) *m){
  /* Calculate statistics for terms have s_functions but not
     others. */
  ETYPE(EXEC_THROUGH_TERMS_INTO)(m, m->workspace, {
      if (mtp->s_func && !mtp->c_func && !mtp->d_func) {
	(*(mtp->s_func))(mtp, nwp);  /* Call d_??? function */
	for (unsigned int k = 0; k < N_CHANGE_STATS; k++) {
	  dstats[k] = mtp->dstats[k]; // Overwrite, not accumulate.
	}
      }
    });
}
