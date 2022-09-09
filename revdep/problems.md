# Blaunet

<details>

* Version: 2.1.0
* GitHub: NA
* Source code: https://github.com/cran/Blaunet
* Date/Publication: 2020-05-22 08:10:11 UTC
* Number of recursive dependencies: 67

Run `revdep_details(, "Blaunet")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘Blaunet-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: active
    > ### Title: Quick summary of blau object.
    > ### Aliases: active
    > 
    > ### ** Examples
    > 
    > data(TwoCities)
    > b <- blau(TwoCities, node.ids = 'respID', ecology.ids = 'samp')
    Error in dimnames(x) <- dn : 
      length of 'dimnames' [1] not equal to array extent
    Calls: blau -> rownames<-
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

* Version: 1.0.1
* GitHub: https://github.com/tidymodels/broom
* Source code: https://github.com/cran/broom
* Date/Publication: 2022-08-29 21:00:08 UTC
* Number of recursive dependencies: 292

Run `revdep_details(, "broom")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'epiR', 'spdep', 'spatialreg'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘spatialreg’, ‘epiR’
    ```

# EpiModel

<details>

* Version: 2.3.0
* GitHub: https://github.com/EpiModel/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2022-07-19 03:30:06 UTC
* Number of recursive dependencies: 106

Run `revdep_details(, "EpiModel")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in ‘EpiModel-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: as.data.frame.netdx
    > ### Title: Extract Timed Edgelists for netdx Objects
    > ### Aliases: as.data.frame.netdx
    > ### Keywords: extract
    > 
    > ### ** Examples
    > 
    ...
    Finished MPLE.
    Stopping at the initial estimate.
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
    +             verbose = FALSE)
    Error in ergm_proposal.NULL(constraints, arguments = if (observational) control$obs.MCMC.prop.args else control$MCMC.prop.args,  : 
      NULL passed to ergm_proposal. This may be due to passing an ergm object from an earlier version. If this is the case, please refit it with the latest version, and try again. If this is not the case, this may be a bug, so please file a bug report.
    Calls: netdx ... eval -> <Anonymous> -> ergm_proposal -> ergm_proposal.NULL
    Execution halted
    ```

*   checking tests ...
    ```
      Running ‘test-all.R’
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
        6.   │   └─base::eval(cl, parent.frame())
        7.   ├─stats::simulate(...)
        8.   └─ergm:::simulate.formula_lhs_network(...)
        9.     ├─ergm::simulate_formula(...)
       10.     ├─tergm:::simulate_formula.network(...)
       11.     │ └─base::eval.parent(mc)
       12.     │   └─base::eval(expr, p)
       13.     │     └─base::eval(expr, p)
       14.     └─ergm (local) `<fn>`(...)
       15.       ├─ergm::ergm_proposal(...)
       16.       └─ergm:::ergm_proposal.NULL(...)
      
      [ FAIL 84 | WARN 55 | SKIP 83 | PASS 589 ]
      Error: Test failures
      Execution halted
    ```

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
      ...
    --- re-building ‘Intro.Rmd’ using rmarkdown
    --- finished re-building ‘Intro.Rmd’
    
    --- re-building ‘attributes-and-summary-statistics.Rmd’ using rmarkdown
    Quitting from lines 307-344 (attributes-and-summary-statistics.Rmd) 
    Error: processing vignette 'attributes-and-summary-statistics.Rmd' failed with diagnostics:
    NULL passed to ergm_proposal. This may be due to passing an ergm object from an earlier version. If this is the case, please refit it with the latest version, and try again. If this is not the case, this may be a bug, so please file a bug report.
    --- failed re-building ‘attributes-and-summary-statistics.Rmd’
    ...
    --- failed re-building ‘model-parameters.Rmd’
    
    --- re-building ‘network-objects.Rmd’ using rmarkdown
    --- finished re-building ‘network-objects.Rmd’
    
    SUMMARY: processing the following files failed:
      ‘attributes-and-summary-statistics.Rmd’ ‘model-parameters.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

# ergMargins

<details>

* Version: 0.1.3
* GitHub: NA
* Source code: https://github.com/cran/ergMargins
* Date/Publication: 2021-06-30 07:40:02 UTC
* Number of recursive dependencies: 58

Run `revdep_details(, "ergMargins")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘methods’ ‘sna’ ‘statnet.common’
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
* Number of recursive dependencies: 65

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
* Number of recursive dependencies: 89

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

# latentnet

<details>

* Version: 2.10.6
* GitHub: https://github.com/statnet/latentnet
* Source code: https://github.com/cran/latentnet
* Date/Publication: 2022-05-11 12:30:05 UTC
* Number of recursive dependencies: 112

Run `revdep_details(, "latentnet")` for more info

</details>

## In both

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# lolog

<details>

* Version: 1.3
* GitHub: https://github.com/statnet/lolog
* Source code: https://github.com/cran/lolog
* Date/Publication: 2021-07-01 07:50:06 UTC
* Number of recursive dependencies: 80

Run `revdep_details(, "lolog")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 24.8Mb
      sub-directories of 1Mb or more:
        libs  23.1Mb
    ```

# netmediate

<details>

* Version: 0.1.0
* GitHub: NA
* Source code: https://github.com/cran/netmediate
* Date/Publication: 2022-08-31 07:50:02 UTC
* Number of recursive dependencies: 99

Run `revdep_details(, "netmediate")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘mediation’
    ```

# sand

<details>

* Version: 2.0.0
* GitHub: https://github.com/kolaczyk/sand
* Source code: https://github.com/cran/sand
* Date/Publication: 2020-07-02 07:20:06 UTC
* Number of recursive dependencies: 160

Run `revdep_details(, "sand")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 6 marked UTF-8 strings
    ```

# stargazer

<details>

* Version: 5.2.3
* GitHub: NA
* Source code: https://github.com/cran/stargazer
* Date/Publication: 2022-03-04 11:50:02 UTC
* Number of recursive dependencies: 0

Run `revdep_details(, "stargazer")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages which this enhances but not available for checking:
      'AER', 'betareg', 'brglm', 'censReg', 'dynlm', 'eha', 'erer',
      'fGarch', 'gee', 'glmx', 'gmm', 'lfe', 'lmtest', 'mclogit', 'mlogit',
      'ordinal', 'plm', 'pscl', 'quantreg', 'rms', 'sampleSelection',
      'spdep'
    ```

# statnetWeb

<details>

* Version: 0.5.6
* GitHub: NA
* Source code: https://github.com/cran/statnetWeb
* Date/Publication: 2020-08-05 18:00:03 UTC
* Number of recursive dependencies: 65

Run `revdep_details(, "statnetWeb")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# texreg

<details>

* Version: 1.38.6
* GitHub: https://github.com/leifeld/texreg
* Source code: https://github.com/cran/texreg
* Date/Publication: 2022-04-06 22:00:02 UTC
* Number of recursive dependencies: 85

Run `revdep_details(, "texreg")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages which this enhances but not available for checking:
      'AER', 'alpaca', 'betareg', 'Bergm', 'bife', 'biglm', 'brglm',
      'brms', 'btergm', 'dynlm', 'eha', 'erer', 'feisr', 'fGarch',
      'fixest', 'forecast', 'gamlss', 'gamlss.inf', 'gee', 'glmmTMB',
      'gmm', 'gnm', 'h2o', 'lfe', 'lqmm', 'maxLik', 'metaSEM', 'mfx',
      'mhurdle', 'miceadds', 'mlogit', 'mnlogit', 'MuMIn', 'oglmx',
      'ordinal', 'pglm', 'plm', 'rms', 'robust', 'simex', 'spatialreg',
      'spdep', 'speedglm', 'truncreg', 'VGAM', 'Zelig'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘h2o’, ‘spatialreg’, ‘eha’, ‘MuMIn’, ‘Bergm’, ‘mfx’, ‘betareg’, ‘bife’, ‘biglm’, ‘brglm’, ‘brms’, ‘btergm’, ‘ordinal’, ‘dynlm’, ‘forecast’, ‘fGarch’, ‘alpaca’, ‘feisr’, ‘lfe’, ‘fixest’, ‘gamlss’, ‘gamlss.inf’, ‘gee’, ‘gmm’, ‘miceadds’, ‘glmmTMB’, ‘gnm’, ‘AER’, ‘robust’, ‘lqmm’, ‘rms’, ‘erer’, ‘maxLik’, ‘mhurdle’, ‘mlogit’, ‘oglmx’, ‘plm’, ‘pglm’, ‘simex’, ‘speedglm’, ‘truncreg’, ‘VGAM’, ‘metaSEM’
    ```

# xergm.common

<details>

* Version: 1.7.8
* GitHub: https://github.com/leifeld/xergm.common
* Source code: https://github.com/cran/xergm.common
* Date/Publication: 2020-04-07 09:50:02 UTC
* Number of recursive dependencies: 35

Run `revdep_details(, "xergm.common")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘RSiena’
    ```

