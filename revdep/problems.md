# broom

<details>

* Version: 0.7.0
* GitHub: https://github.com/tidymodels/broom
* Source code: https://github.com/cran/broom
* Date/Publication: 2020-07-09 12:30:09 UTC
* Number of recursive dependencies: 277

Run `revdep_details(, "broom")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    ...
    NOTE: Versions before 3.6.1 had a bug in the implementation of the bd()
    constraint which distorted the sampled distribution somewhat. In
    addition, Sampson's Monks datasets had mislabeled vertices. See the
    NEWS and the documentation for more details.
    
    NOTE: Some common term arguments pertaining to vertex attribute and
    level selection have changed in 3.10.0. See terms help for more
    details. Use ‘options(ergm.term=list(version="3.9.4"))’ to use old
    behavior.
    
    
    Attaching package: ‘xergm.common’
    
    The following object is masked from ‘package:ergm’:
    
        gof
    
    Loading required package: ggplot2
    Error: package or namespace load failed for ‘btergm’:
     object ‘remove.offset.formula’ is not exported by 'namespace:ergm'
    Execution halted
    ```

# btergm

<details>

* Version: 1.9.9
* GitHub: https://github.com/leifeld/btergm
* Source code: https://github.com/cran/btergm
* Date/Publication: 2020-06-18 05:00:06 UTC
* Number of recursive dependencies: 73

Run `revdep_details(, "btergm")` for more info

</details>

## Newly broken

*   checking whether package ‘btergm’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/homes/morrism/GitHub/StatnetOrganization/ergm/revdep/checks/btergm/new/btergm.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking whether package ‘btergm’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: no DISPLAY variable so Tk is not available
    See ‘/homes/morrism/GitHub/StatnetOrganization/ergm/revdep/checks/btergm/old/btergm.Rcheck/00install.out’ for details.
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
Error: object ‘remove.offset.formula’ is not exported by 'namespace:ergm'
Execution halted
ERROR: lazy loading failed for package ‘btergm’
* removing ‘/homes/morrism/GitHub/StatnetOrganization/ergm/revdep/checks/btergm/new/btergm.Rcheck/btergm’

```
### CRAN

```
* installing *source* package ‘btergm’ ...
** package ‘btergm’ successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Warning: no DISPLAY variable so Tk is not available
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
Warning: no DISPLAY variable so Tk is not available
** testing if installed package can be loaded from final location
Warning: no DISPLAY variable so Tk is not available
** testing if installed package keeps a record of temporary installation path
* DONE (btergm)

```
# EpiModel

<details>

* Version: 2.0.2
* GitHub: https://github.com/statnet/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2020-08-05 09:52:19 UTC
* Number of recursive dependencies: 99

Run `revdep_details(, "EpiModel")` for more info

</details>

## Newly broken

*   checking for missing documentation entries ... WARNING
    ```
    Error: package ‘ergm’ required by ‘tergm’ could not be found
    Call sequence:
    5: stop(gettextf("package %s required by %s could not be found", 
           sQuote(pkg), sQuote(pkgname)), call. = FALSE, domain = NA)
    4: .getRequiredPackages2(pkgInfo, quietly = quietly)
    3: library(pkg, character.only = TRUE, logical.return = TRUE, lib.loc = lib.loc, 
           quietly = quietly)
    2: .getRequiredPackages2(pkgInfo, quietly = quietly)
    1: library(package, lib.loc = lib.loc, character.only = TRUE, verbose = FALSE)
    Execution halted
    All user-level objects in a package should have documentation entries.
    See chapter ‘Writing R documentation files’ in the ‘Writing R
    Extensions’ manual.
    ```

*   checking for code/documentation mismatches ... WARNING
    ```
    ...
    Execution halted
    Error: package ‘ergm’ required by ‘tergm’ could not be found
    Call sequence:
    5: stop(gettextf("package %s required by %s could not be found", 
           sQuote(pkg), sQuote(pkgname)), call. = FALSE, domain = NA)
    4: .getRequiredPackages2(pkgInfo, quietly = quietly)
    3: library(pkg, character.only = TRUE, logical.return = TRUE, lib.loc = lib.loc, 
           quietly = quietly)
    2: .getRequiredPackages2(pkgInfo, quietly = quietly)
    1: library(package, lib.loc = lib.loc, character.only = TRUE, verbose = FALSE)
    Execution halted
    Error: package ‘ergm’ required by ‘tergm’ could not be found
    Call sequence:
    5: stop(gettextf("package %s required by %s could not be found", 
           sQuote(pkg), sQuote(pkgname)), call. = FALSE, domain = NA)
    4: .getRequiredPackages2(pkgInfo, quietly = quietly)
    3: library(pkg, character.only = TRUE, logical.return = TRUE, lib.loc = lib.loc, 
           quietly = quietly)
    2: .getRequiredPackages2(pkgInfo, quietly = quietly)
    1: library(package, lib.loc = lib.loc, character.only = TRUE, verbose = FALSE)
    Execution halted
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘ergm’
    ```

*   checking Rd \usage sections ... NOTE
    ```
    Error: package ‘ergm’ required by ‘tergm’ could not be found
    Call sequence:
    5: stop(gettextf("package %s required by %s could not be found", 
           sQuote(pkg), sQuote(pkgname)), call. = FALSE, domain = NA)
    4: .getRequiredPackages2(pkgInfo, quietly = quietly)
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
    ...
    
    networkDynamic: version 0.10.1, created on 2020-01-16
    Copyright (c) 2020, Carter T. Butts, University of California -- Irvine
                        Ayn Leslie-Cook, University of Washington
                        Pavel N. Krivitsky, University of Wollongong
                        Skye Bender-deMoll, University of Washington
                        with contributions from
                        Zack Almquist, University of California -- Irvine
                        David R. Hunter, Penn State University
                        Li Wang
                        Kirk Li, University of Washington
                        Steven M. Goodreau, University of Washington
                        Jeffrey Horner
                        Martina Morris, University of Washington
    Based on "statnet" project software (statnet.org).
    For license and citation information see statnet.org/attribution
    or type citation("networkDynamic").
    
    Loading required package: tergm
    Error: package ‘ergm’ required by ‘tergm’ could not be found
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
                          Li Wang
                          Kirk Li, University of Washington
                          Steven M. Goodreau, University of Washington
                          Jeffrey Horner
                          Martina Morris, University of Washington
      Based on "statnet" project software (statnet.org).
      For license and citation information see statnet.org/attribution
      or type citation("networkDynamic").
      
      Loading required package: tergm
      Failed with error:  'package 'ergm' required by 'tergm' could not be found'
      Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) : 
        there is no package called 'ergm'
      Calls: test_check ... loadNamespace -> withRestarts -> withOneRestart -> doWithOneRestart
      Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘ndtv’
    ```

# ergm

<details>

* Version: 3.10.4
* GitHub: https://github.com/statnet/ergm
* Source code: https://github.com/cran/ergm
* Date/Publication: 2019-06-10 05:30:07 UTC
* Number of recursive dependencies: 71

Run `revdep_details(, "ergm")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    ...
    Running examples in ‘ergm-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: ergmMPLE
    > ### Title: ERGM Predictors and response for logistic regression calculation
    > ###   of MPLE
    > ### Aliases: ergmMPLE
    > ### Keywords: models regression
    > 
    > ### ** Examples
    > 
    > 
    > data(faux.mesa.high)
    > formula <- faux.mesa.high ~ edges + nodematch("Sex") + nodefactor("Grade")
    > mplesetup <- ergmMPLE(formula)
    Warning: 'compact.rle' is deprecated.
    Use 'compress' instead.
    See help("Deprecated")
    Error in as.rle(x) : could not find function "as.rle"
    Calls: ergmMPLE ... ergm_conlist -> eval -> eval -> <Anonymous> -> rlebdm
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/constrain_degrees_edges.R’ failed.
    Last 13 lines of output:
      > id <- function(nw) apply(as.matrix(nw, matrix.type="adjacency"), 2, sum)
      > e <- function(nw) network.edgecount(nw)
      > 
      > ###### Directed
      > y0 <- as.network(n, density=d, directed=TRUE)
      > 
      > ### Outdegrees
      > ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~odegrees, coef=rep(0,n*2), nsim=nsim, output="stats")
      Error in as.rle(x) : could not find function "as.rle"
      Calls: simulate ... ergm_conlist -> eval -> eval -> <Anonymous> -> rlebdm
      In addition: Warning message:
      'compact.rle' is deprecated.
      Use 'compress' instead.
      See help("Deprecated") 
      Execution halted
    ```

*   checking R code for possible problems ... NOTE
    ```
    rlebdm: no visible global function definition for ‘as.rle’
    Undefined global functions or variables:
      as.rle
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  6.5Mb
      sub-directories of 1Mb or more:
        R      1.2Mb
        doc    1.7Mb
        libs   2.5Mb
    ```

