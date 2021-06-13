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
* Number of recursive dependencies: 85

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
* Number of recursive dependencies: 295

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

*   checking installed package size ... NOTE
    ```
      installed size is 10.5Mb
      sub-directories of 1Mb or more:
        doc    3.6Mb
        help   1.5Mb
        libs   3.9Mb
    ```

# ergm.count

<details>

* Version: 3.4.0
* GitHub: https://github.com/statnet/ergm.count
* Source code: https://github.com/cran/ergm.count
* Date/Publication: 2019-05-15 07:42:59 UTC
* Number of recursive dependencies: 26

Run `revdep_details(, "ergm.count")` for more info

</details>

## Newly broken

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
* Number of recursive dependencies: 59

Run `revdep_details(, "ergmito")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  7.6Mb
      sub-directories of 1Mb or more:
        R      1.1Mb
        libs   5.7Mb
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

## In both

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

*   checking examples ... ERROR
    ```
    Running examples in ‘lolog-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: as.network.Rcpp_DirectedNet
    > ### Title: Convert a DirectedNet to a network object
    > ### Aliases: as.network.Rcpp_DirectedNet
    > 
    > ### ** Examples
    > 
    > el <- matrix(c(1,2),ncol=2)
    > 
    > #make an UndirectedNet with one edge and 5 nodes
    > net <- new(UndirectedNet, el, 5L)
    > 
    > nw <- as.network(net)
    Error in apply(x[, 1:2], 1, sort) : dim(X) must have a positive length
    Calls: as.network ... as.network.matrix -> network.edgelist -> t -> apply
    Execution halted
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 24.9Mb
      sub-directories of 1Mb or more:
        libs  23.1Mb
    ```

# sand

<details>

* Version: 2.0.0
* GitHub: https://github.com/kolaczyk/sand
* Source code: https://github.com/cran/sand
* Date/Publication: 2020-07-02 07:20:06 UTC
* Number of recursive dependencies: 163

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
* Number of recursive dependencies: 53

Run `revdep_details(, "statnetWeb")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# tergmLite

<details>

* Version: 2.2.1
* GitHub: NA
* Source code: https://github.com/cran/tergmLite
* Date/Publication: 2020-07-22 16:50:03 UTC
* Number of recursive dependencies: 67

Run `revdep_details(, "tergmLite")` for more info

</details>

## In both

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
* Number of recursive dependencies: 26

Run `revdep_details(, "xergm.common")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘RSiena’
    ```

