# Bergm

<details>

* Version: 5.0.1
* GitHub: NA
* Source code: https://github.com/cran/Bergm
* Date/Publication: 2019-11-04 12:00:02 UTC
* Number of recursive dependencies: 43

Run `revdep_details(, "Bergm")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    ...
    
    R is free software and comes with ABSOLUTELY NO WARRANTY.
    You are welcome to redistribute it under certain conditions.
    Type 'license()' or 'licence()' for distribution details.
    
      Natural language support but running in an English locale
    
    R is a collaborative project with many contributors.
    Type 'contributors()' for more information and
    'citation()' on how to cite R or R packages in publications.
    
    Type 'demo()' for some demos, 'help()' for on-line help, or
    'help.start()' for an HTML browser interface to help.
    Type 'q()' to quit R.
    
    > pkgname <- "Bergm"
    > source(file.path(R.home("share"), "R", "examples-header.R"))
    > options(warn = 1)
    > library('Bergm')
    Error: package ‘ergm’ required by ‘Bergm’ could not be found
    Execution halted
    ```

## Newly fixed

*   checking examples ... WARNING
    ```
    Found the following significant warnings:
    
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
    Deprecated functions may be defunct as soon as of the next release of
    R.
    See ?Deprecated.
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
    See ‘/home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/btergm/new/btergm.Rcheck/00install.out’ for details.
    ```

## Newly fixed

*   checking whether package ‘btergm’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Warning: no DISPLAY variable so Tk is not available
    See ‘/home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/btergm/old/btergm.Rcheck/00install.out’ for details.
    ```

*   checking examples ... WARNING
    ```
    ...
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
    Deprecated functions may be defunct as soon as of the next release of
    R.
    See ?Deprecated.
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
* removing ‘/home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/btergm/new/btergm.Rcheck/btergm’

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
* Number of recursive dependencies: 100

Run `revdep_details(, "EpiModel")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    ...
    > 
    > # Initialize and parameterize the network model
    > nw <- network_initialize(n = 100)
    > formation <- ~edges
    > target.stats <- 50
    > coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    > 
    > # Model estimation
    > est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    Starting maximum pseudolikelihood estimation (MPLE):
    Evaluating the predictor and response matrix.
    Maximizing the pseudolikelihood.
    Finished MPLE.
    Stopping at the initial estimate.
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE, verbose = FALSE)
    Error in get(name, envir = asNamespace(pkg), inherits = FALSE) : 
      object '.deinf' not found
    Calls: netdx ... stergm_MCMC_sample -> stergm_MCMC_slave -> ::: -> get
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
      ══ testthat results  ═══════════════════════════════════════════════════════════
      [ OK: 242 | SKIPPED: 78 | WARNINGS: 0 | FAILED: 36 ]
      1. Error: Copying attributes from network to attribute list (@test-attr-copy.R#23) 
      2. Error: (unknown) (@test-get.R#18) 
      3. Error: merge for netsim (@test-merge.R#71) 
      4. Error: merge for netsim (@test-merge.R#90) 
      5. Error: merge works for open sims saving nw stats (@test-merge.R#111) 
      6. Error: mutate_epi.netsim (@test-mutate.R#16) 
      7. Error: status.vector and infTime.vector (@test-net-long.R#747) 
      8. Error: tergmLite: 1G, Closed (@test-net-tergmLite.R#22) 
      9. Error: tergmLite: 2G, Closed (@test-net-tergmLite.R#67) 
      1. ...
      
      Error: testthat unit tests failed
      Execution halted
    ```

## Newly fixed

*   checking examples ... WARNING
    ```
    ...
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
    Deprecated functions may be defunct as soon as of the next release of
    R.
    See ?Deprecated.
    ```

# ergm.ego

<details>

* Version: 0.5
* GitHub: https://github.com/statnet/ergm.ego
* Source code: https://github.com/cran/ergm.ego
* Date/Publication: 2019-05-31 16:00:03 UTC
* Number of recursive dependencies: 58

Run `revdep_details(, "ergm.ego")` for more info

</details>

## Newly broken

*   checking Rd cross-references ... WARNING
    ```
    Missing link or links in documentation object 'node-attr-api.Rd':
      ‘node-attr’
    
    See section 'Cross-references' in the 'Writing R Extensions' manual.
    ```

## Newly fixed

*   checking examples ... WARNING
    ```
    ...
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
      Warning: 'compact.rle' is deprecated.
    Deprecated functions may be defunct as soon as of the next release of
    R.
    See ?Deprecated.
    ```

