# Bergm

<details>

* Version: 5.0.2
* GitHub: NA
* Source code: https://github.com/cran/Bergm
* Date/Publication: 2020-11-12 22:20:03 UTC
* Number of recursive dependencies: 37

Run `revdep_details(, "Bergm")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in ‘Bergm-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: bergm
    > ### Title: Parameter estimation for Bayesian ERGMs
    > ### Aliases: bergm
    > 
    > ### ** Examples
    > 
    > # Load the florentine marriage network
    ...
    > # Posterior parameter estimation:
    > p.flo <- bergm(flomarriage ~ edges + kstar(2),
    +                burn.in    = 50,
    +                aux.iters  = 500,
    +                main.iters = 1000,
    +                gamma      = 1.2)
    Error in ergm.Cprepare(y, model) : 
      could not find function "ergm.Cprepare"
    Calls: bergm
    Execution halted
    ```

*   checking R code for possible problems ... NOTE
    ```
    bergm: no visible global function definition for ‘ergm.Cprepare’
    bergmC: no visible global function definition for ‘ergm.Cprepare’
    bergmM: no visible global function definition for ‘ergm.Cprepare’
    ergmAPL: no visible global function definition for ‘ergm.Cprepare’
    Undefined global functions or variables:
      ergm.Cprepare
    ```

# Blaunet

<details>

* Version: 2.1.0
* GitHub: NA
* Source code: https://github.com/cran/Blaunet
* Date/Publication: 2020-05-22 08:10:11 UTC
* Number of recursive dependencies: 84

Run `revdep_details(, "Blaunet")` for more info

</details>

## Newly broken

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

## In both

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

* Version: 0.7.6
* GitHub: https://github.com/tidymodels/broom
* Source code: https://github.com/cran/broom
* Date/Publication: 2021-04-05 20:30:02 UTC
* Number of recursive dependencies: 291

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
* Number of recursive dependencies: 72

Run `revdep_details(, "btergm")` for more info

</details>

## Newly broken

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
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (btergm)


```
# EpiModel

<details>

* Version: 2.0.5
* GitHub: https://github.com/statnet/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2021-05-15 18:50:03 UTC
* Number of recursive dependencies: 102

Run `revdep_details(, "EpiModel")` for more info

</details>

## Newly broken

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

* Version: 3.11.0
* GitHub: https://github.com/statnet/ergm
* Source code: https://github.com/cran/ergm
* Date/Publication: 2020-10-14 09:30:02 UTC
* Number of recursive dependencies: 69

Run `revdep_details(, "ergm")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘Rmpi’
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 10.5Mb
      sub-directories of 1Mb or more:
        doc    3.6Mb
        help   1.5Mb
        libs   3.9Mb
    ```

# ergm.ego

<details>

* Version: 0.6.1
* GitHub: https://github.com/statnet/ergm.ego
* Source code: https://github.com/cran/ergm.ego
* Date/Publication: 2020-11-19 23:10:05 UTC
* Number of recursive dependencies: 55

Run `revdep_details(, "ergm.ego")` for more info

</details>

## Newly broken

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
        Numeric: lengths (784, 792) differ
      Execution halted
    ```

# ergm.rank

<details>

* Version: 1.2.0
* GitHub: https://github.com/statnet/ergm.rank
* Source code: https://github.com/cran/ergm.rank
* Date/Publication: 2019-05-15 07:43:03 UTC
* Number of recursive dependencies: 26

Run `revdep_details(, "ergm.rank")` for more info

</details>

## Newly broken

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
* Number of recursive dependencies: 75

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

# ergmito

<details>

* Version: 0.3-0
* GitHub: https://github.com/muriteams/ergmito
* Source code: https://github.com/cran/ergmito
* Date/Publication: 2020-08-10 21:40:02 UTC
* Number of recursive dependencies: 59

Run `revdep_details(, "ergmito")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in ‘ergmito-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: vcov.ergmito
    > ### Title: Estimation of ERGMs using Maximum Likelihood Estimation (MLE)
    > ### Aliases: vcov.ergmito ergmito
    > 
    > ### ** Examples
    > 
    > 
    ...
    Starting Monte Carlo maximum likelihood estimation (MCMLE):
    Iteration 1 of at most 60:
    Optimizing with step length 1.0000.
    The log-likelihood improved by 0.0007.
    Convergence test p-value: 0.0001. Converged with 99% confidence.
    Finished MCMLE.
    Evaluating log-likelihood at the estimate. Error in UseMethod("%ergmlhs%") : 
      no applicable method for '%ergmlhs%' applied to an object of class "c('matrix', 'array', 'double', 'numeric')"
    Calls: ergm ... ergm.bridge.dindstart.llk -> ergm_preprocess_response -> %ergmlhs%
    Execution halted
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  9.7Mb
      sub-directories of 1Mb or more:
        R      1.1Mb
        libs   7.9Mb
    ```

# gwdegree

<details>

* Version: 0.1.1
* GitHub: https://github.com/michaellevy/gwdegree
* Source code: https://github.com/cran/gwdegree
* Date/Publication: 2016-07-09 10:46:45
* Number of recursive dependencies: 82

Run `revdep_details(, "gwdegree")` for more info

</details>

## Newly broken

*   checking whether package ‘gwdegree’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: replacing previous import ‘ergm::symmetrize’ by ‘sna::symmetrize’ when loading ‘gwdegree’
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/gwdegree/new/gwdegree.Rcheck/00install.out’ for details.
    ```

## In both

*   checking R code for possible problems ... NOTE
    ```
    simCCCent : <anonymous>: no visible global function definition for
      ‘simulate.formula’
    Undefined global functions or variables:
      simulate.formula
    ```

# hergm

<details>

* Version: 4.1-7
* GitHub: NA
* Source code: https://github.com/cran/hergm
* Date/Publication: 2021-03-07 00:00:11 UTC
* Number of recursive dependencies: 66

Run `revdep_details(, "hergm")` for more info

</details>

## Newly broken

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
* Number of recursive dependencies: 75

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
* Number of recursive dependencies: 49

Run `revdep_details(, "NetMix")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  8.3Mb
      sub-directories of 1Mb or more:
        libs   7.9Mb
    ```

# sand

<details>

* Version: 2.0.0
* GitHub: https://github.com/kolaczyk/sand
* Source code: https://github.com/cran/sand
* Date/Publication: 2020-07-02 07:20:06 UTC
* Number of recursive dependencies: 156

Run `revdep_details(, "sand")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 6 marked UTF-8 strings
    ```

# xergm.common

<details>

* Version: 1.7.8
* GitHub: https://github.com/leifeld/xergm.common
* Source code: https://github.com/cran/xergm.common
* Date/Publication: 2020-04-07 09:50:02 UTC
* Number of recursive dependencies: 26

Run `revdep_details(, "xergm.common")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘RSiena’
    ```

