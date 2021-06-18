# Blaunet

<details>

* Version: 2.1.0
* GitHub: NA
* Source code: https://github.com/cran/Blaunet
* Date/Publication: 2020-05-22 08:10:11 UTC
* Number of recursive dependencies: 86

Run `revdep_details(, "Blaunet")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘Blaunet-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: blau
    > ### Title: Converts raw data into an object for Blau status analysis
    > ### Aliases: blau
    > 
    > ### ** Examples
    > 
    > ##simple example
    ...
    > ##example with relational data
    > data(BSANet)
    > square.data <- BSANet$square.data
    > el <- BSANet$el
    > 
    > b <- blau(square.data, node.ids = 'person', ecology.ids = 'city', graph = el)
    Error in Ops.factor(sources, targets) : 
      level sets of factors are different
    Calls: blau ... as.network.data.frame -> .validate_edge_df -> which -> Ops.factor
    Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'gWidgets', 'gWidgetsRGtk2'
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘RGtk2’ ‘cairoDevice’ ‘ergm’ ‘foreign’ ‘haven’ ‘plot3D’ ‘plot3Drgl’
      ‘rgl’ ‘sna’ ‘statnet.common’
      All declared Imports should be used.
    ```

# broom

<details>

* Version: 0.7.7
* GitHub: https://github.com/tidymodels/broom
* Source code: https://github.com/cran/broom
* Date/Publication: 2021-06-13 04:40:17 UTC
* Number of recursive dependencies: 299

Run `revdep_details(, "broom")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'rgeos', 'spdep', 'spatialreg'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘spatialreg’
    ```

# btergm

<details>

* Version: 1.9.13
* GitHub: https://github.com/leifeld/btergm
* Source code: https://github.com/cran/btergm
* Date/Publication: 2020-10-26 14:30:02 UTC
* Number of recursive dependencies: 75

Run `revdep_details(, "btergm")` for more info

</details>

## In both

*   checking whether package ‘btergm’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/btergm/new/btergm.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘btergm’ ...
** package ‘btergm’ successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object ‘ergm.Cprepare’ is not exported by 'namespace:ergm'
Execution halted
ERROR: lazy loading failed for package ‘btergm’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/btergm/new/btergm.Rcheck/btergm’


```
### CRAN

```
* installing *source* package ‘btergm’ ...
** package ‘btergm’ successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error: object ‘ergm.Cprepare’ is not exported by 'namespace:ergm'
Execution halted
ERROR: lazy loading failed for package ‘btergm’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/btergm/old/btergm.Rcheck/btergm’


```
# EpiModel

<details>

* Version: 2.0.5
* GitHub: https://github.com/statnet/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2021-05-15 18:50:03 UTC
* Number of recursive dependencies: 103

Run `revdep_details(, "EpiModel")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘EpiModel-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: as.data.frame.netdx
    > ### Title: Extract Timed Edgelists netdx Objects
    > ### Aliases: as.data.frame.netdx
    > ### Keywords: extract
    > 
    > ### ** Examples
    > 
    ...
    Stopping at the initial estimate.
    Warning: Using x$coef to access the coefficient vector of an ergm is deprecated. Use coef(x) instead.
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
    +             verbose = FALSE)
    Error in is.durational(formation) : 
      could not find function "is.durational"
    Calls: netdx -> simulate -> simulate.network -> is.lasttoggle
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
        3.     └─EpiModel:::FUN(X[[i]], ...)
        4.       └─EpiModel:::netsim_loop(x, param, init, control, s)
        5.         ├─base::withCallingHandlers(...)
        6.         ├─base::do.call(...)
        7.         └─(function (x, param, init, control, s) ...
        8.           └─EpiModel::sim_nets(x, nw, nsteps = control$nsteps, control)
        9.             ├─base::suppressWarnings(...)
       10.             │ └─base::withCallingHandlers(...)
       11.             ├─stats::simulate(...)
       12.             └─tergm:::simulate.network(...)
       13.               └─tergm:::is.lasttoggle(nw, formation, dissolution, monitor)
      
      [ FAIL 38 | WARN 1 | SKIP 79 | PASS 278 ]
      Error: Test failures
      Execution halted
    ```

# ergm

<details>

* Version: 4.0.0
* GitHub: https://github.com/statnet/ergm
* Source code: https://github.com/cran/ergm
* Date/Publication: 2021-06-14 14:30:02 UTC
* Number of recursive dependencies: 81

Run `revdep_details(, "ergm")` for more info

</details>

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      for degeneracy, use the mcmc.diagnostics() function.
      ══ Skipped tests ═══════════════════════════════════════════════════════════════
      • Skipping the rest for time. (1)
      
      ══ Failed tests ════════════════════════════════════════════════════════════════
      ── Failure (test-target-offset.R:23:3): target+offset in a curved ERGM ─────────
      `expect_warning(...)` did not throw the expected warning.
      Backtrace:
          █
       1. └─testthat::expect_warning(...) test-target-offset.R:23:2
       2.   └─testthat:::expect_condition_matching(...)
      
      [ FAIL 1 | WARN 2 | SKIP 1 | PASS 1284 ]
      Error: Test failures
      Execution halted
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  8.3Mb
      sub-directories of 1Mb or more:
        R      1.1Mb
        doc    2.2Mb
        libs   3.8Mb
    ```

# ergm.count

<details>

* Version: 3.4.0
* GitHub: https://github.com/statnet/ergm.count
* Source code: https://github.com/cran/ergm.count
* Date/Publication: 2019-05-15 07:42:59 UTC
* Number of recursive dependencies: 29

Run `revdep_details(, "ergm.count")` for more info

</details>

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/valued_fit.R’ failed.
    Last 13 lines of output:
      
           Null Deviance:   0.00  on 20  degrees of freedom
       Residual Deviance: -17.07  on 19  degrees of freedom
       
      Note that the null model likelihood and deviance are defined to be 0.
      This means that all likelihood-based inference (LRT, Analysis of
      Deviance, AIC, BIC, etc.) is only valid between models with the same
      reference distribution and constraints.
      
      AIC: -15.07  BIC: -14.07  (Smaller is better. MC Std. Err. = 0.3462)
      > true.llk <- sum(dpois(na.omit(c(m)), exp(coef(efit)), log=TRUE)) - sum(dpois(na.omit(c(m)), 1, log=TRUE))
      > 
      > stopifnot(abs(coef(efit)-truth)<0.02)
      Error: abs(coef(efit) - truth) < 0.02 is not TRUE
      Execution halted
    ```

# ergm.ego

<details>

* Version: 0.6.1
* GitHub: https://github.com/statnet/ergm.ego
* Source code: https://github.com/cran/ergm.ego
* Date/Publication: 2020-11-19 23:10:05 UTC
* Number of recursive dependencies: 58

Run `revdep_details(, "ergm.ego")` for more info

</details>

## Newly broken

*   checking Rd \usage sections ... NOTE
    ```
    Error: package ‘ergm’ does not have a namespace and should be re-installed
    Call sequence:
    4: stop(gettextf("package %s does not have a namespace and should be re-installed", 
           sQuote(package)), domain = NA)
    3: library(pkg, character.only = TRUE, logical.return = TRUE, lib.loc = lib.loc, 
           quietly = quietly)
    2: .getRequiredPackages2(pkgInfo, quietly = quietly)
    1: library(package, lib.loc = lib.loc, character.only = TRUE, verbose = FALSE)
    Execution halted
    The \usage entries for S3 methods should use the \method markup and not
    their full name.
    See chapter ‘Writing R documentation files’ in the ‘Writing R
    Extensions’ manual.
    ```

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘ergm.ego-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: simulate.ergm.ego
    > ### Title: Simulate from a 'ergm.ego' fit.
    > ### Aliases: simulate.ergm.ego
    > ### Keywords: models
    > 
    > ### ** Examples
    > 
    ...
    This model was fit using MCMC.  To examine model diagnostics and check
    for degeneracy, use the mcmc.diagnostics() function.
    > colMeans(egosim <- simulate(egofit, popsize=300,nsim=50,
    +                        output="stats", control=control.simulate.ergm.ego(
    +                        simulate.control=control.simulate.formula(MCMC.burnin=2e6))))
    Note: Constructed network has size 205 different from requested 300. Simulated statistics may need to be rescaled.
    Error in simulate.formula_lhs(object = fmh.ego ~ offset(netsize.adj) +  : 
      No applicable method for LHS of type ‘ergm_state_full’, ‘ergm_state_send’, ‘ergm_state_receive’, ‘ergm_state’.
    Calls: colMeans ... eval -> eval -> <Anonymous> -> simulate.formula_lhs
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/EgoStat.tests.R’ failed.
    Last 13 lines of output:
      +   mm("a") + mm("a", levels2=~-1) + mm("a", levels2=-2) + mm("a", levels2=-(2:3)) + mm(~a>7) + mm(a~b) + mm(.~a) + offset(mm(.~a)) + mm("a", levels2 = 1) + 
      +   mm("b", levels = c("a", "c", "e")) + mm("b", levels = c("a", "c", "e"), levels2 = 3) +
      + 
      +   meandeg
      > 
      > f.y <- statnet.common::nonsimp_update.formula(f, y~.)
      > environment(f.y) <- globalenv()
      > f.y.e <- statnet.common::nonsimp_update.formula(f, y.e~.)
      > environment(f.y.e) <- globalenv()
      > 
      > stopifnot(all.equal(summary(f.y),summary(f.y.e)))
      Error: summary(f.y) and summary(f.y.e) are not equal:
        Names: 187 string mismatches
        Numeric: lengths (780, 788) differ
      Execution halted
    ```

# ergm.rank

<details>

* Version: 1.2.0
* GitHub: https://github.com/statnet/ergm.rank
* Source code: https://github.com/cran/ergm.rank
* Date/Publication: 2019-05-15 07:43:03 UTC
* Number of recursive dependencies: 29

Run `revdep_details(, "ergm.rank")` for more info

</details>

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/termTests_rank.R’ failed.
    Last 13 lines of output:
      +                   rank.nonconformity("local1")+
      +                   rank.nonconformity("local2")+
      +                   rank.nonconformity("localAND")+
      +                   rank.deference+
      +                   rank.nodeicov("v")+
      +                   rank.edgecov("m")+
      +                   rank.inconsistency(nw0,"r",xa),
      +                 coef=rep(0,8),response="r", reference=~DiscUnif(1, n-1), nsim=S, statsonly=FALSE)
      Best valid proposal 'DiscUnif' cannot take into account hint(s) 'sparse'.
      Error in as.list(defaultvalues) : 
        argument "response" is missing, with no default
      Calls: simulate ... eval -> eval -> <Anonymous> -> check.ErgmTerm -> as.list
      In addition: Warning message:
      Use of 'statsonly=' argument has been deprecated. Use 'output='stats'' instead. 
      Execution halted
    ```

# ergMargins

<details>

* Version: 0.1.2
* GitHub: NA
* Source code: https://github.com/cran/ergMargins
* Date/Publication: 2021-02-23 14:50:05 UTC
* Number of recursive dependencies: 78

Run `revdep_details(, "ergMargins")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘btergm’ ‘methods’ ‘sna’ ‘statnet’
      All declared Imports should be used.
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘margins’
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# ergmito

<details>

* Version: 0.3-0
* GitHub: https://github.com/muriteams/ergmito
* Source code: https://github.com/cran/ergmito
* Date/Publication: 2020-08-10 21:40:02 UTC
* Number of recursive dependencies: 62

Run `revdep_details(, "ergmito")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  7.1Mb
      sub-directories of 1Mb or more:
        libs   5.7Mb
    ```

# gwdegree

<details>

* Version: 0.1.1
* GitHub: https://github.com/michaellevy/gwdegree
* Source code: https://github.com/cran/gwdegree
* Date/Publication: 2016-07-09 10:46:45
* Number of recursive dependencies: 83

Run `revdep_details(, "gwdegree")` for more info

</details>

## In both

*   checking whether package ‘gwdegree’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: replacing previous import ‘ergm::symmetrize’ by ‘sna::symmetrize’ when loading ‘gwdegree’
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/gwdegree/new/gwdegree.Rcheck/00install.out’ for details.
    ```

*   checking R code for possible problems ... NOTE
    ```
    simCCCent : <anonymous>: no visible global function definition for
      ‘simulate.formula’
    Undefined global functions or variables:
      simulate.formula
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# hergm

<details>

* Version: 4.1-7
* GitHub: NA
* Source code: https://github.com/cran/hergm
* Date/Publication: 2021-03-07 00:00:11 UTC
* Number of recursive dependencies: 70

Run `revdep_details(, "hergm")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    hergm: no visible global function definition for ‘ergm.Cprepare’
    hergm.getnetwork: no visible global function definition for
      ‘ergm.Cprepare’
    hergm.mcmc: no visible global function definition for ‘ergm.Cprepare’
    hergm.preprocess: no visible global function definition for
      ‘ergm.Cprepare’
    hergm.set.mcmc: no visible global function definition for
      ‘ergm.Cprepare’
    Undefined global functions or variables:
      ergm.Cprepare
    ```

# lolog

<details>

* Version: 1.2
* GitHub: https://github.com/statnet/lolog
* Source code: https://github.com/cran/lolog
* Date/Publication: 2019-01-12 22:52:41 UTC
* Number of recursive dependencies: 78

Run `revdep_details(, "lolog")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 33.9Mb
      sub-directories of 1Mb or more:
        libs  32.0Mb
    ```

# NetMix

<details>

* Version: 0.2.0
* GitHub: https://github.com/solivella/NetMix
* Source code: https://github.com/cran/NetMix
* Date/Publication: 2021-03-01 17:40:08 UTC
* Number of recursive dependencies: 52

Run `revdep_details(, "NetMix")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  8.2Mb
      sub-directories of 1Mb or more:
        libs   7.9Mb
    ```

# sand

<details>

* Version: 2.0.0
* GitHub: https://github.com/kolaczyk/sand
* Source code: https://github.com/cran/sand
* Date/Publication: 2020-07-02 07:20:06 UTC
* Number of recursive dependencies: 164

Run `revdep_details(, "sand")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 6 marked UTF-8 strings
    ```

# statnetWeb

<details>

* Version: 0.5.6
* GitHub: NA
* Source code: https://github.com/cran/statnetWeb
* Date/Publication: 2020-08-05 18:00:03 UTC
* Number of recursive dependencies: 54

Run `revdep_details(, "statnetWeb")` for more info

</details>

## In both

*   checking whether package ‘statnetWeb’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: replacing previous import ‘ergm::symmetrize’ by ‘sna::symmetrize’ when loading ‘statnetWeb’
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/statnetWeb/new/statnetWeb.Rcheck/00install.out’ for details.
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# tergm

<details>

* Version: 3.7.0
* GitHub: https://github.com/statnet/tergm
* Source code: https://github.com/cran/tergm
* Date/Publication: 2020-10-15 12:40:03 UTC
* Number of recursive dependencies: 46

Run `revdep_details(, "tergm")` for more info

</details>

## In both

*   checking whether package ‘tergm’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/tergm/new/tergm.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘tergm’ ...
** package ‘tergm’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
gcc -I"/srv/scratch/z3528859/R-4.1.0/include" -DNDEBUG  -I'/srv/scratch/z3528859/R-4.1.0/library/ergm/include' -I/usr/local/include   -fpic  -g -O2  -c DynSA.c -o DynSA.o
In file included from MCMCDyn.h:13,
                 from DynSA.h:13,
                 from DynSA.c:10:
/srv/scratch/z3528859/R-4.1.0/library/ergm/include/edgetree.h:10:9: note: #pragma message: warning: The header file "edgetree.h" has been deprecated in favor of "ergm_edgetree.h" and may be removed in the future.
 #pragma message ("warning: The header file \"edgetree.h\" has been deprecated in favor of \"ergm_edgetree.h\" and may be removed in the future.")
...
                ^~
DynSA.c:92:11: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     if(nwp->duration_info.lasttoggle)
           ^~
DynSA.c:93:27: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     memcpy(lasttoggle, nwp->duration_info.lasttoggle, sizeof(int)*DYADCOUNT(*n_nodes, *bipartite, *dflag));
                           ^~
make: *** [DynSA.o] Error 1
ERROR: compilation failed for package ‘tergm’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/tergm/new/tergm.Rcheck/tergm’


```
### CRAN

```
* installing *source* package ‘tergm’ ...
** package ‘tergm’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
gcc -I"/srv/scratch/z3528859/R-4.1.0/include" -DNDEBUG  -I'/srv/scratch/z3528859/R-4.1.0/library/ergm/include' -I/usr/local/include   -fpic  -g -O2  -c DynSA.c -o DynSA.o
In file included from MCMCDyn.h:13,
                 from DynSA.h:13,
                 from DynSA.c:10:
/srv/scratch/z3528859/R-4.1.0/library/ergm/include/edgetree.h:10:9: note: #pragma message: warning: The header file "edgetree.h" has been deprecated in favor of "ergm_edgetree.h" and may be removed in the future.
 #pragma message ("warning: The header file \"edgetree.h\" has been deprecated in favor of \"ergm_edgetree.h\" and may be removed in the future.")
...
                ^~
DynSA.c:92:11: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     if(nwp->duration_info.lasttoggle)
           ^~
DynSA.c:93:27: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     memcpy(lasttoggle, nwp->duration_info.lasttoggle, sizeof(int)*DYADCOUNT(*n_nodes, *bipartite, *dflag));
                           ^~
make: *** [DynSA.o] Error 1
ERROR: compilation failed for package ‘tergm’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/tergm/old/tergm.Rcheck/tergm’


```
# tergmLite

<details>

* Version: 2.2.1
* GitHub: NA
* Source code: https://github.com/cran/tergmLite
* Date/Publication: 2020-07-22 16:50:03 UTC
* Number of recursive dependencies: 70

Run `revdep_details(, "tergmLite")` for more info

</details>

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      ── Error (test-updateModelTermInputs.R:493:5): gwesp_truedecay ─────────────────
      Error: could not find function "is.durational"
      Backtrace:
          █
       1. └─EpiModel::initialize.net(est, param, init, control, s = 1) test-updateModelTermInputs.R:493:4
       2.   └─EpiModel::sim_nets(x, nw, nsteps = 1, control)
       3.     ├─base::suppressWarnings(...)
       4.     │ └─base::withCallingHandlers(...)
       5.     ├─stats::simulate(...)
       6.     └─tergm:::simulate.network(...)
       7.       └─tergm:::is.lasttoggle(nw, formation, dissolution, monitor)
      
      [ FAIL 24 | WARN 2 | SKIP 0 | PASS 0 ]
      Error: Test failures
      Execution halted
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# xergm.common

<details>

* Version: 1.7.8
* GitHub: https://github.com/leifeld/xergm.common
* Source code: https://github.com/cran/xergm.common
* Date/Publication: 2020-04-07 09:50:02 UTC
* Number of recursive dependencies: 29

Run `revdep_details(, "xergm.common")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘RSiena’
    ```

