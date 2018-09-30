% API for `ergm` Terms
% Pavel N. Krivitsky

# Summary

This document seeks to be the most up-to-date API documentation for `ergm` terms. Note that it is not intended to be a tutorial as much as a description of what inputs and outputs different parts of the system expect.

# Overview

## Types of storage and types of term

The storage API defines two types of storage: *private storage*, which is attached the `ModelTerm` structure and is specific to each `ergm` term, and *public storage*, which is attached to the `Network` and can be accessed by all terms.

A *statistic* is a familiar `ergm` term like "`edges`" or "`nodefactor`": it adds at least one sufficient statistic to the model. Every statistic can have private storage, and it can read from public storage, but it cannot write to public storage.

An *auxiliary* in an `ergm` term but not an ERGM term in the mathematical sense: it adds no statistics to the model and exists only to initialize and maintain public storage to be used by statistics. It may not be specified on an `ergm` formula by the end-user, but only requested by a statistic.

An auxiliary can rely on another auxiliary's public storage. Note that circular dependencies are not checked.

## Code path

For the purposes of this overview, the following information is relevant, and is elaborated formally later:

* Each term has a a pointer `void *mtp->storage` to private storage
* The term has a pointer to an array of pointers to public storage `void **mtp->aux_storage`. (The pointer is to the same location for all terms of a model.)
* Each term has some subset of five functions on the C side: initializers (`i_`), updaters (`u_`), change stats (`c_`), difference stats (`d_`), and finalizers (`f_`). 
* An `InitErgmTerm.` function's output list can have an additional element, `auxiliaries`, a one-sided formula.

1. The user passes a formula to the procedure (`ergm`, `summary`, etc.); this formula may only have statistics.
1. `ergm.getmodel()` is called.
    1. It iterates through the formula terms, calling `InitErgmTerm.<NAME>()` functions in turn, adding their output to the model term list. Some terms include `auxiliaries` formulas in their list.
    1. `update.ErgmTerm()` when it finds that a term has requested auxiliaries, reserves one space for each auxiliary term requested in the requesting term's input vector and shifts `ParamsBeforeCov` accordingly.
    1. It calls `ergm.auxstorage()` with the complete model.
        1. It iterates through initialized model terms, looking at their `auxiliaries` element for a one-sided formula listing their requested auxiliaries.
        1. It iterates through the auxiliary formulas, calling term initialization on them.
        1. It constructs a list of unique initialized auxilariy terms (that is, when two or more terms had requested an identical auxiliary).
        1. It stores the initialized unique auxiliary terms in `model$model.aux` and the index (in the unique list) of the auxiliary requested by each statistic in the reserved space in the input vector of the requesting statistic.
1. `ergm.Cprepare()` is called.
    1. It constructs a `Clist` for the given model
	1. It constructs a `Clist` for `model.aux`.
	1. It concatenates the two `Clist`s. (The number of terms includes both statistics and auxiliaries.)
1. The result is passed to the C code.
1. `ModelInitialize()` initializes all terms (statistic and auxiliaries), also counting up the number of auxiliaries (distinguished by having no `c_`, `d_`, or `s_` functions). If a term exports a `c_` function, a `d_` function is not looked up.
1. `ModelInitialize()` counts up the number of auxiliaries and allocates an array of `void *`s, one for each auxiliary. The `mtp->aux_storage` pointer for each term is set to point to that (one) array.
1. `NetworkInitialize()` initializes the network.
1. `InitStats()` is called, calling the initializer (`i_` function) of each term or, if not found, an updater (`u_` function) with invalid input (i.e., toggle $(0,0)$) is called in case the term developer prefers a one-function implementation.
    * The terms are initialized in reverse order, so auxiliaries are initialized before the statistics are, and statistics can count on their auxiliaries being initialized by the time they are initialized.
1. For each iteration, a proposal is made.
    1. `ChangeStats()` is called.
        1. It calls the `d_` functions, for those terms for which they are initialized.
	    1. It iterates through the proposed toggles.
	        1. It calls the `c_` functions.
		    1. It adds its output to the cumulative change.
		    1. It calls the `u_` functions with the toggle (if more to come).
		    1. It makes the toggle provisionally (if more to come).
	    1. It undoes the provisional toggles.
    1. If the proposal is accepted, `u_` function is called and network is updated for each toggle.
1. `DestroyStats()` is called, iterating through the terms.
    1. Finalizer (`f_`) function is called, if defined.
	1. If `mtp->storage` is not `NULL`, it is freed.
1. `ModelDestroy()` is called, freeing elements of `mtp->aux_storage` if not `NULL`.
1. `NetworkDestroy()` is called.

# Formal API definition

## New struct elements

`void *mtp->storage`: A pointer to private storage.

`void **mtp->aux_storage`: An array of pointers to public storage, referring to the same memory location for all terms in the model.

# Evaluation API

## `R` side

Unchanged.

## `C` side

### `c_`, `d_`, and `s_` functions

`d_` functions are the original difference statistics. `c_` functions are new, while `s_` functions have been around for a long time, but never formally documented.

Change statistic: `void c_<NAME>(Vertex tail, Vertex head, Model *mtp, Network *nwp)`

Difference statistic: `void d_<NAME>(Vertex *tails, Vertex *heads, Model *mtp, Network *nwp)`

Summary statistic: `void s_<NAME>(Model *mtp, Network *nwp)`

#### Expect
##### Parameters

`Edge ntgoggles`: Number of edges to be toggled.

`Vertex tail`: Tail of (1) dyad to be toggled.

`Vertex *tails`: An array of tails of the dyads to be toggled.

`Vertex head`: Head of (1) dyad to be toggled.

`Vertex *heads`: An array of heads of the dyads to be toggled.

`Model *mtp`: A pointer to a `Model` of interest.

`Network *nwp`: A pointer to the `Network` of interest before any toggles are applied.

##### Storage

All functions except for `s_` expect any storage they need to be initialized and up to date (consistent with `nwp`). In particular, if their statistic requested $k$ auxiliary terms, the first $k$ elements of its `mtp->inputparam` (`INPUT_PARAM`) vector will be the indexes of `mtp->aux_storage` where they can find the respective objects, with the rest of their inputs shifted over.

#### Macros

It is worth noting that macros defined for `d_` functions that refer to a specific toggle, such as `TAIL`, `HEAD`, etc. might not be usable in a `c_` function, but it's made up for by `c_` function's reduced need for bookkeeping: `tail`, `head`, etc. can be used directly.

#### Side-effects

These functions overwrite `mtp->dstats` (often aliased as `CHANGE_STAT`) with the following:

* `c_` and `d_` functions: change of the value of the statistic they implement relative to `nwp` due to the toggles.
* `s_` the value of the statistic it implements.

# Storage API

Every `ergm` term has private storage, found at `void *mtp->storage`, which allows it to store arbitrary information about the state of the network, as well as precalculated values of variables, preallocated memory it needs for its calculations, or any other use. It does so by specifying an updating function (and, optionally, an initialization and a finalization function). This updating function is called every time the network is about to change. The API for these functions is defined below.

Public storage is found at `void **mtp->aux_storage`. Each auxiliary term gets assigned a slot (i.e., `void *nwp->mtp->aux_storage[i]`) to manage; its slot number is the first element of its input vector, and terms requesting it are told which slot to look in in a similar fashion. An auxiliary term that requests other auxiliaries will have its own slot as the first input and the slots of auxiliaries it requests as subsequent inputs.

## `R` side

A statistic that only references its private storage *or* is an auxiliary itself does not need to do anything special on the `R` side.

To request an auxiliary, a term's `InitErgmTerm` call's output list must include an `auxiliaries` element containing a one-sided `ergm`-style formula listing the auxiliary terms it wishes to use separated by the `+` operator.

## `C` side: Modifying Storage

### `i_` functions: Initializer/Constructor

This function is optional for using storage: if it's not provided, the model code will call the `u_` function with an invalid toggle first, signaling for it to initialize.

`void i_<NAME>(Model *mtp, Network *nwp)`

#### Expects

In general, `i_` function expects to be called after `ModelInitialize()` and `NetworkInitialize()`, before any `c_` or `d_` functions. That is, the network must be populated with the ties of its initial state and have `mtp->aux_storage` vector allocated.

##### Private storage

Network populated with initial ties and initialized model.

##### Public storage

The first element of `mtp->inputparams` (a.k.a. `INPUT_PARAM`) is the index of the element of `mtp->aux_storage` to be managed by this auxiliary. That is `mtp->aux_storage[(int)INPUT_PARAM[0]]` is a `void *` to point to the data to be public.

The other data passed from the `InitErgmTerm.` are shifted over to make room for it.

#### Side-effects
##### Private storage
Allocates memory for the information to be stored and overwrites `mtp->storage` with a pointer to it, then updates the stored information to be consistent with `*nwp`.

###### Public storage

Allocates memory for the information to be stored and overwrites `mtp->aux_storage[(int)INPUT_PARAM[0]]` with a pointer to it, then updates the stored information to be consistent with `*nwp`.

An auxiliary can also use its private storage as needed.

### `u_` functions: Updater

`void u_<NAME>(Vertex tail, Vertex head, Model *mtp, Network *nwp)`

#### Expects

Initialized network. If no `i_` function was provided, to be called with a $(0,0)$ toggle as a signal to initialize; otherwise, initialized storage. Any statistic can rely on auxiliaries being initialized before it.

#### Side-effects

If called with an invalid toggle and uninitialized storage, initialize.

Update the state of its storage (`mtp->storage` and/or `mtp->aux_storage[(int)INPUT_PARAM[0]]`) to match what the state of `*nwp` would be after the given dyad had been toggled.

### `f_` functions: Finalizer/Destructor

This function is optional for using storage: if it's not provided, the model code will free any pointers to `mtp->aux_storage` and `mtp->storage` that are not `NULL`.

`void f_<NAME>(Model *mtp, Network *nwp)`

#### Expects:

Network and a model.

#### Side-effects

Deallocates its storage (`mtp->storage` and/or `mtp->aux_storage[(int)INPUT_PARAM[0]]`) and sets its pointers to `NULL`.

## `C` side: Accessing Storage

### Private storage

`c_`, `d_`, and `s_` functions can read from, but not write to, their private storage. `c_` and `d_` functions can rely on initialization having been called before.

### Public storage

Auxilaries must not implement `c_`, `d_`, and `s_` functions.

Terms requesting one or more auxiliaries will be passed the indices of the element of `mtp->aux_storage` by inserting them at the start of `mtp->inputparams` (a.k.a. `INPUT_PARAM`). That is `mtp->aux_storage[(int)INPUT_PARAM[0]]` is a `void *` to point to the data public by the first auxiliary term on the `auxiliaries` formula, `mtp->aux_storage[(int)INPUT_PARAM[1]]` is the second, etc..

# Macros

The following helper macros have been defined to date, and can be found in `storage.h`.

## Memory management

`ALLOC_STORAGE(nmemb, stored_type, store_into)`: Allocate a vector of `nmemb` elements of type `stored_type`, save its pointer to private storage and also to a `stored_type *store_into` which is also declared. Should be used by the `i_` function, but may also be used by the `u_` function.

`GET_STORAGE(stored_type, store_into)`: Declare `stored_type *store_into` and assign the pointer to private storage to it. Can be used by all functions.


`ALLOC_AUX_STORAGE(nmemb, stored_type, store_into)`: Allocate a vector of `nmemb` elements of type `stored_type`, and save it to the auxiliary storage slot belonging to the calling auxiliary and into a `stored_type *store_into` which is also declared. Can be used by the `i_` function, but may also be used by the `u_` function.

`GET_AUX_STORAGE(stored_type, store_into)`: Declare `stored_type *store_into` and assign the pointer to the auxiliary storage (either for a statistic or for the auxiliary). Can bn used by all functions.

`GET_AUX_STORAGE_NUM(stored_type, store_into, ind)`: Declare `stored_type *store_into` and assign the pointer to the `ind`th auxiliary). Can be used by clients of auxiliaries.

## Miscellaneous helpers

`ALLOC_AUX_SOCIOMATRIX(stored_type, store_into)`: Allocate an array of appropriate dimension with elements of type `stored_type`, save it to auxiliary storage, and into `**store_type`, so that `store_into[i][j]` returns the value associated with dyad $(i,j)$, with vertices indexed from 1. For bipartite and undirected networks, as little space as possible (resp. a rectangle or a triangle) is allocated.

Note that this term assumes that the private and the public storage of the calling term are not used in any other way.

`FREE_AUX_SOCIOMATRIX`: Frees the sociomatrix allocated by `ALLOC_AUX_SOCIOMATRIX`.

# MHproposal storage API

## Auxiliaries

MHproposals may also request auxiliary terms. An `InitMHp.` with an `auxiliaries` formula will similarly receive the positions of its auxiliaries' slots in the network. However, this appears to have a slight cost in speed and a potentially significant cost in memory, since the auxiliary may need to duplicate the information in the `MH_` function.

## Private storage

Functions prefixed with `Mi_`, `Mu_`, and, `Mf_` serve as respectively the initializers, the updaters, and the finalizers of the MHproposal storage, though the old-style call with `MHp->ntoggles==0` is also supported. Macros in the `MHstorage.h` header file can be used to access storage the same way as for the statistics.

The function called to generate the proposal can have a prefix of either `MH_` (for backwards compatibility) or `Mp_` for consistency.

One important difference is that `Mp_` function *is* permitted to write to its private storage. This may be useful if, say, a systematic sample is desired.
