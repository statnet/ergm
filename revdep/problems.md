# Blaunet

<details>

* Version: 2.1.0
* GitHub: NA
* Source code: https://github.com/cran/Blaunet
* Date/Publication: 2020-05-22 08:10:11 UTC
* Number of recursive dependencies: 70

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

* Version: 0.7.8
* GitHub: https://github.com/tidymodels/broom
* Source code: https://github.com/cran/broom
* Date/Publication: 2021-06-24 08:50:02 UTC
* Number of recursive dependencies: 302

Run `revdep_details(, "broom")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘broom-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: glance.Rchoice
    > ### Title: Glance at a(n) Rchoice object
    > ### Aliases: glance.Rchoice
    > 
    > ### ** Examples
    > 
    > library(Rchoice)
    ...
      term         estimate std.error statistic p.value
      <chr>           <dbl>     <dbl>     <dbl>   <dbl>
    1 constant       6.01      4.64     1.30      0.195
    2 mpg           -0.0957    0.134   -0.716     0.474
    3 hp            -0.0244    0.0214  -1.14      0.253
    4 factor(cyl)6  -0.980     1.23    -0.797     0.426
    5 factor(cyl)8  -6.39    938.      -0.00682   0.995
    > glance(mod)
    Error: 'getSummary.Rchoice' is not an exported object from 'namespace:Rchoice'
    Execution halted
    ```

*   checking dependencies in R code ... WARNING
    ```
    Missing or unexported object: ‘Rchoice::getSummary.Rchoice’
    ```

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'rgeos', 'spdep', 'spatialreg'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘spatialreg’
    ```

# ergm

<details>

* Version: 4.0.1
* GitHub: https://github.com/statnet/ergm
* Source code: https://github.com/cran/ergm
* Date/Publication: 2021-06-21 07:20:02 UTC
* Number of recursive dependencies: 81

Run `revdep_details(, "ergm")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  8.3Mb
      sub-directories of 1Mb or more:
        R      1.1Mb
        doc    2.2Mb
        libs   3.8Mb
    ```

# ergMargins

<details>

* Version: 0.1.3
* GitHub: NA
* Source code: https://github.com/cran/ergMargins
* Date/Publication: 2021-06-30 07:40:02 UTC
* Number of recursive dependencies: 56

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

* Version: 1.3
* GitHub: https://github.com/statnet/lolog
* Source code: https://github.com/cran/lolog
* Date/Publication: 2021-07-01 07:50:06 UTC
* Number of recursive dependencies: 78

Run `revdep_details(, "lolog")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 24.8Mb
      sub-directories of 1Mb or more:
        libs  23.1Mb
    ```

# sand

<details>

* Version: 2.0.0
* GitHub: https://github.com/kolaczyk/sand
* Source code: https://github.com/cran/sand
* Date/Publication: 2020-07-02 07:20:06 UTC
* Number of recursive dependencies: 166

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

