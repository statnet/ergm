/*  File src/model.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include <string.h>
#include "model.h"

/*****************
  void ModelDestroy
******************/
void ModelDestroy(Model *m, Network *nwp)
{  
  DestroyStats(nwp, m);
  
  for(unsigned int i=0; i < m->n_aux; i++)
    if(m->termarray[0].aux_storage[i]!=NULL){
      free(m->termarray[0].aux_storage[i]);
      m->termarray[0].aux_storage[i] = NULL;
    }
  
  if(m->termarray[0].aux_storage!=NULL){
    free(m->termarray[0].aux_storage);
  }
  
  EXEC_THROUGH_TERMS({
      if(mtp->aux_storage!=NULL)
	mtp->aux_storage=NULL;
    });
  
  free(m->dstatarray);
  free(m->termarray);
  free(m->workspace); 
  free(m);
}

/*****************
 int ModelInitialize

 Allocate and initialize the ModelTerm structures, each of which contains
 all necessary information about how to compute one term in the model.
*****************/
Model* ModelInitialize (char *fnames, char *sonames, double **inputsp,
			int n_terms) {
  int i, j, k, l, offset;
  ModelTerm *thisterm;
  char *fn,*sn;
  Model *m;
  double *inputs=*inputsp;
  
  m = (Model *) malloc(sizeof(Model));
  m->n_terms = n_terms;
  m->termarray = (ModelTerm *) malloc(sizeof(ModelTerm) * n_terms);
  m->dstatarray = (double **) malloc(sizeof(double *) * n_terms);
  m->n_stats = 0;
  m->n_aux = 0;
  for (l=0; l < n_terms; l++) {
      thisterm = m->termarray + l;
      
      /* Initialize storage and term functions to NULL. */
      thisterm->storage = NULL;
      thisterm->aux_storage = NULL;
      thisterm->d_func = NULL;
      thisterm->c_func = NULL;
      thisterm->s_func = NULL;
      thisterm->i_func = NULL;
      thisterm->u_func = NULL;
      thisterm->f_func = NULL;
      
      /* First, obtain the term name and library: fnames points to a
      single character string, consisting of the names of the selected
      options concatenated together and separated by spaces.  This is
      passed by the calling R function.  These names are matched with
      their respective C functions that calculate the appropriate
      statistics.  Similarly, sonames points to a character string
      containing the names of the shared object files associated with
      the respective functions.*/
      for (; *fnames == ' ' || *fnames == 0; fnames++);
      for (i = 0; fnames[i] != ' ' && fnames[i] != 0; i++);
      fnames[i] = 0;
      for (; *sonames == ' ' || *sonames == 0; sonames++);
      for (j = 0; sonames[j] != ' ' && sonames[j] != 0; j++);
      sonames[j] = 0;
      /* Extract the required string information from the relevant sources */
      if((fn=(char *)malloc(sizeof(char)*(i+3)))==NULL){
        error("Error in ModelInitialize: Can't allocate %d bytes for fn. Memory has not been deallocated, so restart R sometime soon.\n",
		sizeof(char)*(i+3));
      }
      fn[1]='_';
      for(k=0;k<i;k++)
        fn[k+2]=fnames[k];
      fn[i+2]='\0';
      /* fn is now the string 'd_[name]', where [name] is fname */
/*      Rprintf("fn: %s\n",fn); */
      if((sn=(char *)malloc(sizeof(char)*(j+1)))==NULL){
        error("Error in ModelInitialize: Can't allocate %d bytes for sn. Memory has not been deallocated, so restart R sometime soon.\n",
		sizeof(char)*(j+1));
      }
      sn=strncpy(sn,sonames,j);
      sn[j]='\0';


      /*  Second, process the values in
          model$option[[optionnumber]]$inputs; See comments in
          InitErgm.r for details. This needs to be done before change
          statistica are found, to determine whether a term is
          auxiliary.  */
      offset = (int) *inputs++;  /* Set offset for attr vector */
      /*      Rprintf("offsets: %f %f %f %f %f\n",inputs[0],inputs[1],inputs[2], */
      /*		         inputs[3],inputs[4],inputs[5]); */
      thisterm->nstats = (int) *inputs++; /* If >0, # of statistics returned. If ==0 an auxiliary statistic. */
      
      /*      Rprintf("l %d offset %d thisterm %d\n",l,offset,thisterm->nstats); */
      
      /*  Update the running total of statistics */
      m->n_stats += thisterm->nstats; 
      m->dstatarray[l] = (double *) malloc(sizeof(double) * thisterm->nstats);
      thisterm->dstats = m->dstatarray[l];  /* This line is important for
                                               eventually freeing up malloc'ed
					       memory, since thisterm->dstats
					       can be modified but 
					       m->dstatarray[l] cannot be.  */
      thisterm->statcache = (double *) malloc(sizeof(double) * thisterm->nstats);

      thisterm->ninputparams = (int) *inputs++; /* Set # of inputs */
      /* thisterm->inputparams is a ptr to inputs */
      thisterm->inputparams = (thisterm->ninputparams ==0) ? 0 : inputs; 
      
      thisterm->attrib = inputs + offset; /* Ptr to attributes */
      inputs += thisterm->ninputparams;  /* Skip to next model option */

      /* If the term's nstats==0, it is auxiliary: it does not affect
	 acceptance probabilities or contribute any
	 statistics. Therefore, its d_ and s_ functions are never
	 called and are not initialized. It is only used for its u_
	 function. Therefore, the following code is only run when
	 thisterm->nstats>0. */
      if(thisterm->nstats){
	/*  Most important part of the ModelTerm:  A pointer to a
	    function that will compute the change in the network statistic of 
	    interest for a particular edge toggle.  This function is obtained by
	    searching for symbols associated with the object file with prefix
	    sn, having the name fn.  Assuming that one is found, we're golden.*/ 
	fn[0]='c';
	thisterm->c_func = 
	  (void (*)(Vertex, Vertex, ModelTerm*, Network*))
	  R_FindSymbol(fn,sn,NULL);

	if(thisterm->c_func==NULL){
	  fn[0]='d';
	  thisterm->d_func = 
	    (void (*)(Edge, Vertex*, Vertex*, ModelTerm*, Network*))
	    R_FindSymbol(fn,sn,NULL);
	
	  if(thisterm->d_func==NULL){
	    error("Error in ModelInitialize: could not find function %s in "
		  "namespace for package %s. Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
	  }
	}
	
	/* Optional function to compute the statistic of interest for
	   the network given. It can be more efficient than going one
	   edge at a time. */
	fn[0]='s';
	thisterm->s_func = 
	  (void (*)(ModelTerm*, Network*)) R_FindSymbol(fn,sn,NULL);	
      }else m->n_aux++;
      
      /* Optional functions to store persistent information about the
	 network state between calls to d_ functions. */
      fn[0]='u';
      thisterm->u_func = 
	(void (*)(Vertex, Vertex, ModelTerm*, Network*)) R_FindSymbol(fn,sn,NULL);

      /* If it's an auxiliary, then it needs a u_function, or
	 it's not doing anything. */
      if(thisterm->nstats==0 && thisterm->u_func==NULL){
	error("Error in ModelInitialize: could not find updater function %s in "
	      "namespace for package %s: this term will not do anything. Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
      }
  

      /* Optional-optional functions to initialize and finalize the
	 term's storage. */
      
      fn[0]='i';
      thisterm->i_func = 
	(void (*)(ModelTerm*, Network*)) R_FindSymbol(fn,sn,NULL);

      fn[0]='f';
      thisterm->f_func = 
	(void (*)(ModelTerm*, Network*)) R_FindSymbol(fn,sn,NULL);

      
      /*Clean up by freeing sn and fn*/
      free((void *)fn);
      free((void *)sn);

      /*  The lines above set thisterm->inputparams to point to needed input
      parameters (or zero if none) and then increments the inputs pointer so
      that it points to the inputs for the next model option for the next pass
      through the loop. */

      fnames += i;
      sonames += j;
    }
  
  m->workspace = (double *) malloc(sizeof(double) * m->n_stats);
  for(i=0; i < m->n_stats; i++)
    m->workspace[i] = 0.0;

  
  /* Allocate auxiliary storage and put a pointer to it on every model term. */
  if(m->n_aux){
    m->termarray[0].aux_storage = (void *) calloc(m->n_aux, sizeof(void *));
    for(l=1; l < n_terms; l++)
      m->termarray[l].aux_storage = m->termarray[0].aux_storage;
  }
  
  *inputsp = inputs;
  return m;
}

/*
  ChangeStats
  A helper's helper function to compute change statistics.
  The vector of changes is written to m->workspace.
*/
void ChangeStats(unsigned int ntoggles, Vertex *tails, Vertex *heads,
				 Network *nwp, Model *m){
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */ 

  /* Make a pass through terms with d_functions. */
  EXEC_THROUGH_TERMS_INTO(m->workspace, {
      mtp->dstats = dstats; /* Stuck the change statistic here.*/
      if(mtp->c_func==NULL && mtp->d_func)
	(*(mtp->d_func))(ntoggles, tails, heads, 
			 mtp, nwp);  /* Call d_??? function */
    });
  /* Notice that mtp->dstats now points to the appropriate location in
     m->workspace. */
  
  /* Put the original destination dstats back unless there's only one
     toggle. */
  if(ntoggles!=1){
    unsigned int i = 0;
    EXEC_THROUGH_TERMS({
	mtp->dstats = m->dstatarray[i];
	i++;
      });
  }

  /* Make a pass through terms with c_functions. */
  int toggle;
  FOR_EACH_TOGGLE(toggle){
    EXEC_THROUGH_TERMS_INTO(m->workspace, {
	if(mtp->c_func){
	  (*(mtp->c_func))(*(tails+toggle), *(heads+toggle),
			   mtp, nwp);  /* Call d_??? function */
	  
	  if(ntoggles!=1){
	    for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	      dstats[k] += mtp->dstats[k];
	    }
	  }
	}
      });

    /* Execute storage updates */
    IF_MORE_TO_COME(toggle){
      UPDATE_STORAGE_COND(tails[toggle],heads[toggle], m, nwp, mtp->d_func==NULL);
      TOGGLE(tails[toggle],heads[toggle]);
    }
  }
  /* Undo previous storage updates and toggles */
  UNDO_PREVIOUS(toggle){
    UPDATE_STORAGE_COND(tails[toggle],heads[toggle], m, nwp, mtp->d_func==NULL);
    TOGGLE(tails[toggle],heads[toggle]);
  }
}
      
/*
  InitStats
  A helper's helper function to initialize storage for functions that use it.
*/
void InitStats(Network *nwp, Model *m){
  // Iterate in reverse, so that auxliary terms get initialized first.  
  EXEC_THROUGH_TERMS_INREVERSE({
      double *dstats = mtp->dstats;
      mtp->dstats = NULL; // Trigger segfault if i_func tries to write to change statistics.
      if(mtp->i_func)
	(*(mtp->i_func))(mtp, nwp);  /* Call i_??? function */
      else if(mtp->u_func) /* No initializer but an updater -> uses a 1-function implementation. */
	(*(mtp->u_func))(0, 0, mtp, nwp);  /* Call u_??? function */
      mtp->dstats = dstats;
    });
}

/*
  DestroyStats
  A helper's helper function to finalize storage for functions that use it.
*/
void DestroyStats(Network *nwp, Model *m){
  unsigned int i=0;
  EXEC_THROUGH_TERMS({
      if(mtp->f_func)
	(*(mtp->f_func))(mtp, nwp);  /* Call f_??? function */
      free(m->dstatarray[i]);
      free(mtp->statcache);
      if(mtp->storage){
	free(mtp->storage);
	mtp->storage = NULL;
      }
      i++;
    });
}

