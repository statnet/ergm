# Blaunet

<details>

* Version: 2.2.1
* GitHub: NA
* Source code: https://github.com/cran/Blaunet
* Date/Publication: 2022-09-27 08:10:08 UTC
* Number of recursive dependencies: 75

Run `revdepcheck::revdep_details(, "Blaunet")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘digest’ ‘ergm’ ‘foreign’ ‘gWidgets2’ ‘gWidgets2tcltk’ ‘haven’
      ‘plot3D’ ‘plot3Drgl’ ‘rgl’ ‘sna’ ‘statnet.common’
      All declared Imports should be used.
    ```

# broom

<details>

* Version: 1.0.5
* GitHub: https://github.com/tidymodels/broom
* Source code: https://github.com/cran/broom
* Date/Publication: 2023-06-09 22:50:02 UTC
* Number of recursive dependencies: 305

Run `revdepcheck::revdep_details(, "broom")` for more info

</details>

## In both

*   checking tests ...
    ```
      Running ‘spelling.R’
      Running ‘test-all.R’
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
           ▆
        1. └─lme4::lmer(mpg ~ wt + (1 | cyl), data = mtcars) at test-stats-anova.R:133:3
        2.   ├─base::do.call(...)
        3.   └─lme4 (local) `<fn>`(...)
        4.     ├─base::do.call(...)
    ...
        5.     └─methods (local) `<rfMthdDf>`(...)
        6.       └─methods::new(def, ...)
        7.         ├─methods::initialize(value, ...)
        8.         └─methods::initialize(value, ...)
        9.           └─.Object$initialize(...)
       10.             └─lme4 (local) initializePtr()
      
      [ FAIL 1 | WARN 0 | SKIP 81 | PASS 1027 ]
      Error: Test failures
      Execution halted
    ```

# btergm

<details>

* Version: 1.10.11
* GitHub: https://github.com/leifeld/btergm
* Source code: https://github.com/cran/btergm
* Date/Publication: 2023-10-05 20:10:02 UTC
* Number of recursive dependencies: 100

Run `revdepcheck::revdep_details(, "btergm")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 89 marked UTF-8 strings
    ```

# dnr

<details>

* Version: 0.3.5
* GitHub: NA
* Source code: https://github.com/cran/dnr
* Date/Publication: 2020-11-30 17:10:03 UTC
* Number of recursive dependencies: 73

Run `revdepcheck::revdep_details(, "dnr")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

*   checking Rd files ... NOTE
    ```
    checkRd: (-1) dnr.Rd:12: Escaped LaTeX specials: \#
    ```

# ergmito

<details>

* Version: 0.3-1
* GitHub: https://github.com/muriteams/ergmito
* Source code: https://github.com/cran/ergmito
* Date/Publication: 2023-06-14 10:42:05 UTC
* Number of recursive dependencies: 69

Run `revdepcheck::revdep_details(, "ergmito")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  7.3Mb
      sub-directories of 1Mb or more:
        libs   5.9Mb
    ```

# latentnet

<details>

* Version: 2.10.6
* GitHub: https://github.com/statnet/latentnet
* Source code: https://github.com/cran/latentnet
* Date/Publication: 2022-05-11 12:30:05 UTC
* Number of recursive dependencies: 116

Run `revdepcheck::revdep_details(, "latentnet")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘ergm.userterms’
    ```

*   checking Rd files ... NOTE
    ```
    checkRd: (-1) simulate.ergmm.Rd:35: Escaped LaTeX specials: \$
    checkRd: (-1) simulate.ergmm.Rd:36: Escaped LaTeX specials: \$
    ```

*   checking Rd cross-references ... NOTE
    ```
    Unknown package ‘ergm.userterms’ in Rd xrefs
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# lolog

<details>

* Version: 1.3.1
* GitHub: https://github.com/statnet/lolog
* Source code: https://github.com/cran/lolog
* Date/Publication: 2023-12-07 12:40:02 UTC
* Number of recursive dependencies: 86

Run `revdepcheck::revdep_details(, "lolog")` for more info

</details>

## In both

*   R CMD check timed out
    

*   checking installed package size ... NOTE
    ```
      installed size is 27.8Mb
      sub-directories of 1Mb or more:
        libs  26.2Mb
    ```

# sand

<details>

* Version: 2.0.0
* GitHub: https://github.com/kolaczyk/sand
* Source code: https://github.com/cran/sand
* Date/Publication: 2020-07-02 07:20:06 UTC
* Number of recursive dependencies: 156

Run `revdepcheck::revdep_details(, "sand")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘networkTomography’
    ```

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

Run `revdepcheck::revdep_details(, "stargazer")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages which this enhances but not available for checking:
      'brglm', 'censReg', 'dynlm', 'eha', 'erer', 'fGarch', 'glmx',
      'mclogit', 'rms', 'sampleSelection'
    ```

# statnet

<details>

* Version: 2019.6
* GitHub: https://github.com/statnet/statnet
* Source code: https://github.com/cran/statnet
* Date/Publication: 2019-06-14 08:00:06 UTC
* Number of recursive dependencies: 99

Run `revdepcheck::revdep_details(, "statnet")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘networksis’
    ```

*   checking startup messages can be suppressed ... NOTE
    ```
    unable to reach CRAN
    
    It looks like this package (or a package it requires) has a startup
    message which cannot be suppressed: see ?packageStartupMessage.
    ```

# statnetWeb

<details>

* Version: 0.5.6
* GitHub: NA
* Source code: https://github.com/cran/statnetWeb
* Date/Publication: 2020-08-05 18:00:03 UTC
* Number of recursive dependencies: 67

Run `revdepcheck::revdep_details(, "statnetWeb")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# texreg

<details>

* Version: 1.39.3
* GitHub: https://github.com/leifeld/texreg
* Source code: https://github.com/cran/texreg
* Date/Publication: 2023-11-10 00:00:03 UTC
* Number of recursive dependencies: 104

Run `revdepcheck::revdep_details(, "texreg")` for more info

</details>

## In both

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
           ▆
        1. └─lme4::lmer(...) at test-texreg.R:169:3
        2.   ├─base::do.call(...)
        3.   └─lme4 (local) `<fn>`(...)
        4.     ├─base::do.call(...)
        5.     └─methods (local) `<rfMthdDf>`(...)
        6.       └─methods::new(def, ...)
        7.         ├─methods::initialize(value, ...)
        8.         └─methods::initialize(value, ...)
        9.           └─.Object$initialize(...)
       10.             └─lme4 (local) initializePtr()
      
      [ FAIL 1 | WARN 0 | SKIP 32 | PASS 213 ]
      Error: Test failures
      Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Packages which this enhances but not available for checking:
      'alpaca', 'bife', 'brglm', 'brms', 'dynlm', 'eha', 'erer', 'feisr',
      'fGarch', 'forecast', 'gamlss', 'gamlss.inf', 'glmmTMB', 'gnm',
      'h2o', 'logitr', 'lqmm', 'metaSEM', 'mhurdle', 'miceadds', 'mnlogit',
      'MuMIn', 'oglmx', 'pglm', 'rms', 'simex', 'truncreg', 'Zelig'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘h2o’, ‘eha’, ‘MuMIn’, ‘bife’, ‘brglm’, ‘brms’, ‘dynlm’, ‘forecast’, ‘fGarch’, ‘alpaca’, ‘feisr’, ‘gamlss’, ‘gamlss.inf’, ‘miceadds’, ‘glmmTMB’, ‘gnm’, ‘lqmm’, ‘rms’, ‘erer’, ‘mhurdle’, ‘oglmx’, ‘pglm’, ‘simex’, ‘truncreg’, ‘metaSEM’
    ```

