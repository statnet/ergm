# Amelia

<details>

* Version: 1.8.1
* GitHub: NA
* Source code: https://github.com/cran/Amelia
* Date/Publication: 2022-11-19 19:40:29 UTC
* Number of recursive dependencies: 48

Run `revdep_details(, "Amelia")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  5.9Mb
      sub-directories of 1Mb or more:
        doc    1.4Mb
        libs   3.8Mb
    ```

# amt

<details>

* Version: 0.2.1.0
* GitHub: https://github.com/jmsigner/amt
* Source code: https://github.com/cran/amt
* Date/Publication: 2023-03-28 08:20:02 UTC
* Number of recursive dependencies: 149

Run `revdep_details(, "amt")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘amt-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: time_of_day
    > ### Title: Time of the day when a fix was taken
    > ### Aliases: time_of_day time_of_day.track_xyt time_of_day.steps_xyt
    > 
    > ### ** Examples
    > 
    > data(deer)
    > deer |> time_of_day()
    Error in C_force_tz(to_posixct(time), tz, roll_dst) : 
      CCTZ: Unrecognized output timezone: ":/etc/localtime"
    Calls: time_of_day ... <Anonymous> -> .force_tz -> from_posixct -> C_force_tz
    Execution halted
    ```

# AovBay

<details>

* Version: 0.1.0
* GitHub: NA
* Source code: https://github.com/cran/AovBay
* Date/Publication: 2021-07-22 06:30:02 UTC
* Number of recursive dependencies: 144

Run `revdep_details(, "AovBay")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 41.3Mb
      sub-directories of 1Mb or more:
        libs  41.0Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘RcppParallel’
      All declared Imports should be used.
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# archetyper

<details>

* Version: 0.1.0
* GitHub: https://github.com/mkorvink/archetyper
* Source code: https://github.com/cran/archetyper
* Date/Publication: 2021-03-17 14:00:05 UTC
* Number of recursive dependencies: 181

Run `revdep_details(, "archetyper")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘bannerCommenter’ ‘config’ ‘feather’ ‘here’ ‘knitr’ ‘log4r’ ‘ps’
      ‘rmarkdown’ ‘skimr’ ‘snakecase’ ‘testthat’ ‘tidyverse’
      All declared Imports should be used.
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# autocogs

<details>

* Version: 0.1.4
* GitHub: https://github.com/schloerke/autocogs
* Source code: https://github.com/cran/autocogs
* Date/Publication: 2021-05-29 17:00:05 UTC
* Number of recursive dependencies: 74

Run `revdep_details(, "autocogs")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘broom’ ‘diptest’ ‘ggplot2’ ‘hexbin’ ‘MASS’ ‘moments’
      All declared Imports should be used.
    ```

# BasketballAnalyzeR

<details>

* Version: 0.5.0
* GitHub: https://github.com/sndmrc/BasketballAnalyzeR
* Source code: https://github.com/cran/BasketballAnalyzeR
* Date/Publication: 2020-06-26 09:00:11 UTC
* Number of recursive dependencies: 74

Run `revdep_details(, "BasketballAnalyzeR")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘circlize’ ‘hexbin’ ‘scales’ ‘sna’
      All declared Imports should be used.
    ```

# BayesPostEst

<details>

* Version: 0.3.2
* GitHub: https://github.com/ShanaScogin/BayesPostEst
* Source code: https://github.com/cran/BayesPostEst
* Date/Publication: 2021-11-11 08:10:05 UTC
* Number of recursive dependencies: 159

Run `revdep_details(, "BayesPostEst")` for more info

</details>

## In both

*   checking package dependencies ... ERROR
    ```
    Packages required but not available: 'R2jags', 'runjags', 'rjags'
    
    See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
    manual.
    ```

# Blaunet

<details>

* Version: 2.2.1
* GitHub: NA
* Source code: https://github.com/cran/Blaunet
* Date/Publication: 2022-09-27 08:10:08 UTC
* Number of recursive dependencies: 75

Run `revdep_details(, "Blaunet")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘digest’ ‘ergm’ ‘foreign’ ‘gWidgets2’ ‘gWidgets2tcltk’ ‘haven’
      ‘plot3D’ ‘plot3Drgl’ ‘rgl’ ‘sna’ ‘statnet.common’
      All declared Imports should be used.
    ```

# broom.mixed

<details>

* Version: 0.2.9.4
* GitHub: https://github.com/bbolker/broom.mixed
* Source code: https://github.com/cran/broom.mixed
* Date/Publication: 2022-04-17 17:42:29 UTC
* Number of recursive dependencies: 164

Run `revdep_details(, "broom.mixed")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking: 'glmmADMB', 'R2jags'
    ```

# btergm

<details>

* Version: 1.10.10
* GitHub: https://github.com/leifeld/btergm
* Source code: https://github.com/cran/btergm
* Date/Publication: 2023-04-24 12:40:14 UTC
* Number of recursive dependencies: 100

Run `revdep_details(, "btergm")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      
      ══ Failed tests ════════════════════════════════════════════════════════════════
      ── Error ('test-btergm.R:73:3'): offset argument in btergm works without composition change ──
      Error in `ergm::ergm.pl(nw, Clist.miss, model, theta.offset = c(rep(FALSE, 
          length(l$rhs.terms) - 1), TRUE), verbose = FALSE, control = control.ergm)`: VECTOR_ELT() can only be applied to a 'list', not a 'NULL'
      Backtrace:
          ▆
       1. ├─base::suppressWarnings(...) at test-btergm.R:73:2
       2. │ └─base::withCallingHandlers(...)
       3. └─btergm::btergm(...)
       4.   └─ergm::ergm.pl(...)
      
      [ FAIL 1 | WARN 0 | SKIP 8 | PASS 51 ]
      Error: Test failures
      Execution halted
    ```

# card

<details>

* Version: 0.1.0
* GitHub: https://github.com/asshah4/card
* Source code: https://github.com/cran/card
* Date/Publication: 2020-09-03 07:52:10 UTC
* Number of recursive dependencies: 147

Run `revdep_details(, "card")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘Hmisc’ ‘methods’ ‘readr’ ‘recipes’ ‘sf’ ‘stringr’ ‘utils’
      All declared Imports should be used.
    ```

# catfun

<details>

* Version: 0.1.4
* GitHub: NA
* Source code: https://github.com/cran/catfun
* Date/Publication: 2019-06-14 14:10:03 UTC
* Number of recursive dependencies: 112

Run `revdep_details(, "catfun")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘broom’ ‘epitools’
      All declared Imports should be used.
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# convergEU

<details>

* Version: 0.5.4
* GitHub: NA
* Source code: https://github.com/cran/convergEU
* Date/Publication: 2023-02-18 08:50:02 UTC
* Number of recursive dependencies: 175

Run `revdep_details(, "convergEU")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 1 marked UTF-8 string
    ```

# csdata

<details>

* Version: 2023.5.22
* GitHub: https://github.com/csids/csdata
* Source code: https://github.com/cran/csdata
* Date/Publication: 2023-05-22 10:10:18 UTC
* Number of recursive dependencies: 136

Run `revdep_details(, "csdata")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘geojsonio’
    ```

# describedata

<details>

* Version: 0.1.0
* GitHub: NA
* Source code: https://github.com/cran/describedata
* Date/Publication: 2019-08-02 11:50:02 UTC
* Number of recursive dependencies: 69

Run `revdep_details(, "describedata")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# didimputation

<details>

* Version: 0.3.0
* GitHub: NA
* Source code: https://github.com/cran/didimputation
* Date/Publication: 2022-08-25 20:02:33 UTC
* Number of recursive dependencies: 64

Run `revdep_details(, "didimputation")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

# disto

<details>

* Version: 0.2.0
* GitHub: https://github.com/talegari/disto
* Source code: https://github.com/cran/disto
* Date/Publication: 2018-08-02 12:50:02 UTC
* Number of recursive dependencies: 124

Run `revdep_details(, "disto")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘dplyr’ ‘proxy’
      All declared Imports should be used.
    ```

# dnr

<details>

* Version: 0.3.5
* GitHub: NA
* Source code: https://github.com/cran/dnr
* Date/Publication: 2020-11-30 17:10:03 UTC
* Number of recursive dependencies: 73

Run `revdep_details(, "dnr")` for more info

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

# dotwhisker

<details>

* Version: 0.7.4
* GitHub: https://github.com/fsolt/dotwhisker
* Source code: https://github.com/cran/dotwhisker
* Date/Publication: 2021-09-02 14:50:35 UTC
* Number of recursive dependencies: 74

Run `revdep_details(, "dotwhisker")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Unknown package ‘broomExtra’ in Rd xrefs
    ```

# dplyr

<details>

* Version: 1.1.2
* GitHub: https://github.com/tidyverse/dplyr
* Source code: https://github.com/cran/dplyr
* Date/Publication: 2023-04-20 14:00:03 UTC
* Number of recursive dependencies: 96

Run `revdep_details(, "dplyr")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking: 'RMySQL', 'RPostgreSQL'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘RMySQL’
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 4 marked UTF-8 strings
    ```

# dragon

<details>

* Version: 1.2.1
* GitHub: https://github.com/sjspielman/dragon
* Source code: https://github.com/cran/dragon
* Date/Publication: 2022-04-08 08:42:33 UTC
* Number of recursive dependencies: 130

Run `revdep_details(, "dragon")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘htmltools’
      All declared Imports should be used.
    ```

# eoffice

<details>

* Version: 0.2.2
* GitHub: NA
* Source code: https://github.com/cran/eoffice
* Date/Publication: 2022-10-05 07:30:02 UTC
* Number of recursive dependencies: 107

Run `revdep_details(, "eoffice")` for more info

</details>

## In both

*   checking whether package ‘eoffice’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/eoffice/new/eoffice.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘eoffice’ ...
** package ‘eoffice’ successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/srv/scratch/z3528859/R-4.3.0/library/magick/libs/magick.so':
  /apps/z_install_tree/linux-rocky8-ivybridge/gcc-12.2.0/glib-2.74.1-jrm3i3hfvjbyqoo4hipwrurprh3uteor/lib/libgobject-2.0.so.0: undefined symbol: g_uri_ref
Calls: <Anonymous> ... asNamespace -> loadNamespace -> library.dynam -> dyn.load
Execution halted
ERROR: lazy loading failed for package ‘eoffice’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/eoffice/new/eoffice.Rcheck/eoffice’


```
### CRAN

```
* installing *source* package ‘eoffice’ ...
** package ‘eoffice’ successfully unpacked and MD5 sums checked
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/srv/scratch/z3528859/R-4.3.0/library/magick/libs/magick.so':
  /apps/z_install_tree/linux-rocky8-ivybridge/gcc-12.2.0/glib-2.74.1-jrm3i3hfvjbyqoo4hipwrurprh3uteor/lib/libgobject-2.0.so.0: undefined symbol: g_uri_ref
Calls: <Anonymous> ... asNamespace -> loadNamespace -> library.dynam -> dyn.load
Execution halted
ERROR: lazy loading failed for package ‘eoffice’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/eoffice/old/eoffice.Rcheck/eoffice’


```
# EpiModel

<details>

* Version: 2.3.2
* GitHub: https://github.com/EpiModel/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2023-02-16 23:30:02 UTC
* Number of recursive dependencies: 120

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
    Maximizing the pseudolikelihood.
    Finished MPLE.
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
    +             verbose = FALSE)
    Error in simulate.formula_lhs(object = TARGET_STATS ~ edges, nsim = if (dynamic ==  : 
      No applicable method for LHS of type ‘NULL’.
    Calls: netdx ... eval -> eval -> <Anonymous> -> simulate.formula_lhs
    Execution halted
    ```

*   checking tests ...
    ```
      Running ‘test-all.R’
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
       16.             ├─stats::simulate(...)
       17.             └─ergm:::simulate.formula_lhs(...)
      ── Error ('test-utils.R:38:3'): color_tea ──────────────────────────────────────
      Error in `UseMethod("get.vertex.attribute")`: no applicable method for 'get.vertex.attribute' applied to an object of class "NULL"
      Backtrace:
          ▆
       1. └─EpiModel::netsim(est, param, init, control) at test-utils.R:38:2
       2.   └─EpiModel::crosscheck.net(x, param, init, control)
       3.     ├─get_vertex_attribute(nw, "status") %in% c("s", "i", "r")
       4.     └─EpiModel::get_vertex_attribute(nw, "status")
       5.       └─network::get.vertex.attribute(...)
      
      [ FAIL 37 | WARN 37 | SKIP 119 | PASS 376 ]
      Error: Test failures
      Execution halted
    ```

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
    --- re-building ‘attributes-and-summary-statistics.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    
    Quitting from lines 303-338 [tracker-list] (attributes-and-summary-statistics.Rmd)
    Error: processing vignette 'attributes-and-summary-statistics.Rmd' failed with diagnostics:
    no applicable method for 'get.vertex.attribute' applied to an object of class "NULL"
    --- failed re-building ‘attributes-and-summary-statistics.Rmd’
    
    --- re-building ‘Intro.Rmd’ using rmarkdown
    ...
    <p>The cumulative edgelist refers to the historical list of edges in a network with the time step they start and stopped. Such a list allows to query current relationships (contacts, partnerships, etc.) as well as past ones.</p>
    <h3 id="using-the-cumulative-edgelist">Using the Cumulative Edgelist</h3>
    <p>The creation and update of the cumulative edgelist is done through the <code>EpiModel [... truncated]
    --- finished re-building ‘network-objects.Rmd’
    
    SUMMARY: processing the following files failed:
      ‘attributes-and-summary-statistics.Rmd’ ‘model-parameters.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

# ERSA

<details>

* Version: 0.1.3
* GitHub: NA
* Source code: https://github.com/cran/ERSA
* Date/Publication: 2020-09-22 23:00:02 UTC
* Number of recursive dependencies: 99

Run `revdep_details(, "ERSA")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# estimatr

<details>

* Version: 1.0.0
* GitHub: https://github.com/DeclareDesign/estimatr
* Source code: https://github.com/cran/estimatr
* Date/Publication: 2022-07-04 12:40:05 UTC
* Number of recursive dependencies: 99

Run `revdep_details(, "estimatr")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 29.6Mb
      sub-directories of 1Mb or more:
        libs  29.2Mb
    ```

# eurostat

<details>

* Version: 3.8.2
* GitHub: https://github.com/rOpenGov/eurostat
* Source code: https://github.com/cran/eurostat
* Date/Publication: 2023-03-06 15:40:02 UTC
* Number of recursive dependencies: 102

Run `revdep_details(, "eurostat")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 1054 marked UTF-8 strings
    ```

# ExpDes.pt

<details>

* Version: 1.2.2
* GitHub: NA
* Source code: https://github.com/cran/ExpDes.pt
* Date/Publication: 2021-10-05 04:20:08 UTC
* Number of recursive dependencies: 1

Run `revdep_details(, "ExpDes.pt")` for more info

</details>

## In both

*   checking Rd files ... NOTE
    ```
    checkRd: (-1) ccF.Rd:5: Escaped LaTeX specials: \&
    checkRd: (-1) ccF.Rd:32: Escaped LaTeX specials: \&
    checkRd: (-1) ccF.Rd:29: Escaped LaTeX specials: \&
    ```

# export

<details>

* Version: 0.3.0
* GitHub: https://github.com/tomwenseleers/export
* Source code: https://github.com/cran/export
* Date/Publication: 2022-12-07 15:40:02 UTC
* Number of recursive dependencies: 97

Run `revdep_details(, "export")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘export-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: graph2bitmap
    > ### Title: Save currently active R graph to bitmap format
    > ### Aliases: graph2bitmap graph2png graph2tif graph2jpg
    > 
    > ### ** Examples
    > 
    > # Create a file name
    ...
    > graph2tif(x = x, file = filen, dpi = 400, height = 5, aspectr = 4)
    TIFFOpen: /scratch/pbs.4420393.kman.restech.unsw.edu.tif: Permission denied.
    Warning in dev.off(ii) :
      unable to open TIFF file '/scratch/pbs.4420393.kman.restech.unsw.edu.tif'
    Exported graph as /scratch/pbs.4420393.kman.restech.unsw.edu.tif
    > graph2jpg(x = x, file = filen, dpi = 400, height = 5, aspectr = 4)
    Error in grid.newpage() : 
      could not open file '/scratch/pbs.4420393.kman.restech.unsw.edu.jpeg'
    Calls: graph2jpg ... graph2bitmap -> myplot -> print -> print.ggplot -> grid.newpage
    Execution halted
    ```

# eyetrackingR

<details>

* Version: 0.2.0
* GitHub: https://github.com/samhforbes/eyetrackingR
* Source code: https://github.com/cran/eyetrackingR
* Date/Publication: 2021-09-27 10:00:14 UTC
* Number of recursive dependencies: 97

Run `revdep_details(, "eyetrackingR")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘rlang’
      All declared Imports should be used.
    ```

# finalfit

<details>

* Version: 1.0.6
* GitHub: https://github.com/ewenharrison/finalfit
* Source code: https://github.com/cran/finalfit
* Date/Publication: 2023-01-14 13:40:02 UTC
* Number of recursive dependencies: 157

Run `revdep_details(, "finalfit")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  5.6Mb
      sub-directories of 1Mb or more:
        doc   4.9Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘tidyselect’
      All declared Imports should be used.
    ```

# fivethirtyeight

<details>

* Version: 0.6.2
* GitHub: https://github.com/rudeboybert/fivethirtyeight
* Source code: https://github.com/cran/fivethirtyeight
* Date/Publication: 2021-10-07 13:40:02 UTC
* Number of recursive dependencies: 72

Run `revdep_details(, "fivethirtyeight")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘fivethirtyeightdata’
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  5.0Mb
      sub-directories of 1Mb or more:
        data   3.8Mb
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘fivethirtyeightdata’
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 584 marked UTF-8 strings
    ```

# flextable

<details>

* Version: 0.9.1
* GitHub: https://github.com/davidgohel/flextable
* Source code: https://github.com/cran/flextable
* Date/Publication: 2023-04-02 11:20:02 UTC
* Number of recursive dependencies: 138

Run `revdep_details(, "flextable")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘flextable-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: save_as_html
    > ### Title: Save flextable objects in an 'HTML' file
    > ### Aliases: save_as_html
    > 
    > ### ** Examples
    > 
    > ft1 <- flextable(head(iris))
    > tf1 <- tempfile(fileext = ".html")
    > save_as_html(ft1, path = tf1)
    Error in tmp_rmd(title = paste0(title, collapse = "\n"), lang = lang) : 
      pandoc_available() is not TRUE
    Calls: save_as_html -> render_htmltag -> tmp_rmd -> stopifnot
    Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'equatags', 'doconv', 'pdftools'
    ```

# forestmodel

<details>

* Version: 0.6.2
* GitHub: NA
* Source code: https://github.com/cran/forestmodel
* Date/Publication: 2020-07-19 11:50:03 UTC
* Number of recursive dependencies: 58

Run `revdep_details(, "forestmodel")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# fwildclusterboot

<details>

* Version: 0.13.0
* GitHub: https://github.com/s3alfisc/fwildclusterboot
* Source code: https://github.com/cran/fwildclusterboot
* Date/Publication: 2023-02-26 01:00:13 UTC
* Number of recursive dependencies: 124

Run `revdep_details(, "fwildclusterboot")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  6.5Mb
      sub-directories of 1Mb or more:
        libs   5.4Mb
    ```

# ggformula

<details>

* Version: 0.10.4
* GitHub: https://github.com/ProjectMOSAIC/ggformula
* Source code: https://github.com/cran/ggformula
* Date/Publication: 2023-04-11 14:50:02 UTC
* Number of recursive dependencies: 133

Run `revdep_details(, "ggformula")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘akima’
    ```

# ghypernet

<details>

* Version: 1.1.0
* GitHub: NA
* Source code: https://github.com/cran/ghypernet
* Date/Publication: 2021-10-15 13:30:05 UTC
* Number of recursive dependencies: 98

Run `revdep_details(, "ghypernet")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘methods’
      All declared Imports should be used.
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 11 marked UTF-8 strings
    ```

# glmmfields

<details>

* Version: 0.1.7
* GitHub: https://github.com/seananderson/glmmfields
* Source code: https://github.com/cran/glmmfields
* Date/Publication: 2023-03-10 22:50:02 UTC
* Number of recursive dependencies: 109

Run `revdep_details(, "glmmfields")` for more info

</details>

## Newly broken

*   checking installed package size ... NOTE
    ```
      installed size is 75.2Mb
      sub-directories of 1Mb or more:
        libs  74.4Mb
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

## Newly fixed

*   checking whether package ‘glmmfields’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/glmmfields/old/glmmfields.Rcheck/00install.out’ for details.
    ```

# glmmTMB

<details>

* Version: 1.1.7
* GitHub: https://github.com/glmmTMB/glmmTMB
* Source code: https://github.com/cran/glmmTMB
* Date/Publication: 2023-04-05 21:03:31 UTC
* Number of recursive dependencies: 160

Run `revdep_details(, "glmmTMB")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 72.9Mb
      sub-directories of 1Mb or more:
        doc             1.3Mb
        libs           67.7Mb
        test_data       2.2Mb
        vignette_data   1.1Mb
    ```

# goldfish

<details>

* Version: 1.6.4
* GitHub: https://github.com/snlab-ch/goldfish
* Source code: https://github.com/cran/goldfish
* Date/Publication: 2022-08-24 08:30:06 UTC
* Number of recursive dependencies: 153

Run `revdep_details(, "goldfish")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 11.0Mb
      sub-directories of 1Mb or more:
        libs   9.7Mb
    ```

# greybox

<details>

* Version: 1.0.8
* GitHub: https://github.com/config-i1/greybox
* Source code: https://github.com/cran/greybox
* Date/Publication: 2023-04-02 15:30:07 UTC
* Number of recursive dependencies: 72

Run `revdep_details(, "greybox")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package which this enhances but not available for checking: ‘vars’
    ```

# gWQS

<details>

* Version: 3.0.4
* GitHub: NA
* Source code: https://github.com/cran/gWQS
* Date/Publication: 2021-05-20 09:30:02 UTC
* Number of recursive dependencies: 103

Run `revdep_details(, "gWQS")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘dplyr’
      All declared Imports should be used.
    ```

# hbal

<details>

* Version: 1.2.8
* GitHub: NA
* Source code: https://github.com/cran/hbal
* Date/Publication: 2023-01-24 15:30:05 UTC
* Number of recursive dependencies: 86

Run `revdep_details(, "hbal")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  8.6Mb
      sub-directories of 1Mb or more:
        libs   8.3Mb
    ```

# healthyR

<details>

* Version: 0.2.1
* GitHub: https://github.com/spsanderson/healthyR
* Source code: https://github.com/cran/healthyR
* Date/Publication: 2023-04-06 22:20:03 UTC
* Number of recursive dependencies: 154

Run `revdep_details(, "healthyR")` for more info

</details>

## In both

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
      ...
    --- re-building ‘getting-started.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    
    Quitting from lines 87-94 [alos_plot_interactive] (getting-started.Rmd)
    Error: processing vignette 'getting-started.Rmd' failed with diagnostics:
    invalid 'path' argument
    --- failed re-building ‘getting-started.Rmd’
    
    SUMMARY: processing the following file failed:
      ‘getting-started.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  5.5Mb
      sub-directories of 1Mb or more:
        data   1.5Mb
        doc    3.8Mb
    ```

# heplots

<details>

* Version: 1.4-2
* GitHub: https://github.com/friendly/heplots
* Source code: https://github.com/cran/heplots
* Date/Publication: 2022-10-19 23:35:06 UTC
* Number of recursive dependencies: 112

Run `revdep_details(, "heplots")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘Sleuth2’, ‘archdata’, ‘qqtest’
    ```

# highcharter

<details>

* Version: 0.9.4
* GitHub: https://github.com/jbkunst/highcharter
* Source code: https://github.com/cran/highcharter
* Date/Publication: 2022-01-03 16:40:05 UTC
* Number of recursive dependencies: 150

Run `revdep_details(, "highcharter")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘geojsonio’
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 11 marked UTF-8 strings
    ```

# huxtable

<details>

* Version: 5.5.2
* GitHub: https://github.com/hughjonesd/huxtable
* Source code: https://github.com/cran/huxtable
* Date/Publication: 2022-12-16 13:30:02 UTC
* Number of recursive dependencies: 171

Run `revdep_details(, "huxtable")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘R6’ ‘xml2’
      All declared Imports should be used.
    ```

# IceSat2R

<details>

* Version: 1.0.4
* GitHub: https://github.com/mlampros/IceSat2R
* Source code: https://github.com/cran/IceSat2R
* Date/Publication: 2022-11-17 07:50:02 UTC
* Number of recursive dependencies: 144

Run `revdep_details(, "IceSat2R")` for more info

</details>

## In both

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
    --- re-building ‘IceSat-2_Atlas_products_PDF.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    
    Quitting from lines 605-608 [unnamed-chunk-20] (IceSat-2_Atlas_products_PDF.Rmd)
    Error: processing vignette 'IceSat-2_Atlas_products_PDF.Rmd' failed with diagnostics:
    invalid 'path' argument
    --- failed re-building ‘IceSat-2_Atlas_products_PDF.Rmd’
    
    --- re-building ‘IceSat-2_Mission_Orbits_PDF.Rmd’ using rmarkdown
    ...
    <div class="figure" style="text-align: center">
    <img src="himalayas_aoi.png" alt="Display the area of interest in Himalayas" width="105%" height="100%" />
    <p class="capt [... truncated]
    --- finished re-building ‘IceSat-2_Virtual_File_System_Orbits_PDF.Rmd’
    
    SUMMARY: processing the following files failed:
      ‘IceSat-2_Atlas_products_PDF.Rmd’ ‘IceSat-2_Mission_Orbits_PDF.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 10 marked UTF-8 strings
    ```

# industRial

<details>

* Version: 0.1.0
* GitHub: https://github.com/J-Ramalho/industRial
* Source code: https://github.com/cran/industRial
* Date/Publication: 2021-06-11 09:40:02 UTC
* Number of recursive dependencies: 197

Run `revdep_details(, "industRial")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 2 marked UTF-8 strings
    ```

# interactions

<details>

* Version: 1.1.5
* GitHub: https://github.com/jacob-long/interactions
* Source code: https://github.com/cran/interactions
* Date/Publication: 2021-07-02 07:00:04 UTC
* Number of recursive dependencies: 99

Run `revdep_details(, "interactions")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘pequod’
    ```

# ivx

<details>

* Version: 1.1.0
* GitHub: https://github.com/kvasilopoulos/ivx
* Source code: https://github.com/cran/ivx
* Date/Publication: 2020-11-24 13:50:02 UTC
* Number of recursive dependencies: 78

Run `revdep_details(, "ivx")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

# konfound

<details>

* Version: 0.4.0
* GitHub: https://github.com/jrosen48/konfound
* Source code: https://github.com/cran/konfound
* Date/Publication: 2021-06-01 07:40:05 UTC
* Number of recursive dependencies: 141

Run `revdep_details(, "konfound")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘mice’ ‘tibble’
      All declared Imports should be used.
    ```

# latentnet

<details>

* Version: 2.10.6
* GitHub: https://github.com/statnet/latentnet
* Source code: https://github.com/cran/latentnet
* Date/Publication: 2022-05-11 12:30:05 UTC
* Number of recursive dependencies: 115

Run `revdep_details(, "latentnet")` for more info

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

# lcsm

<details>

* Version: 0.3.2
* GitHub: https://github.com/milanwiedemann/lcsm
* Source code: https://github.com/cran/lcsm
* Date/Publication: 2023-02-25 23:40:02 UTC
* Number of recursive dependencies: 136

Run `revdep_details(, "lcsm")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘cli’
      All declared Imports should be used.
    ```

# LDAShiny

<details>

* Version: 0.9.3
* GitHub: https://github.com/JavierDeLaHoz/LDAShiny
* Source code: https://github.com/cran/LDAShiny
* Date/Publication: 2021-03-29 10:02:12 UTC
* Number of recursive dependencies: 139

Run `revdep_details(, "LDAShiny")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘beepr’ ‘broom’ ‘chinese.misc’ ‘dplyr’ ‘DT’ ‘highcharter’
      ‘htmlwidgets’ ‘ldatuning’ ‘parallel’ ‘plotly’ ‘purrr’ ‘quanteda’
      ‘shinyalert’ ‘shinyBS’ ‘shinycssloaders’ ‘shinyjs’ ‘shinyWidgets’
      ‘SnowballC’ ‘stringr’ ‘textmineR’ ‘tidyr’ ‘tidytext’ ‘topicmodels’
      All declared Imports should be used.
    ```

# lin.eval

<details>

* Version: 0.1.2
* GitHub: NA
* Source code: https://github.com/cran/lin.eval
* Date/Publication: 2019-02-22 00:00:03 UTC
* Number of recursive dependencies: 29

Run `revdep_details(, "lin.eval")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# lolog

<details>

* Version: 1.3
* GitHub: https://github.com/statnet/lolog
* Source code: https://github.com/cran/lolog
* Date/Publication: 2021-07-01 07:50:06 UTC
* Number of recursive dependencies: 85

Run `revdep_details(, "lolog")` for more info

</details>

## Newly fixed

*   R CMD check timed out
    

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 27.8Mb
      sub-directories of 1Mb or more:
        libs  26.1Mb
    ```

# lspline

<details>

* Version: 1.0-0
* GitHub: NA
* Source code: https://github.com/cran/lspline
* Date/Publication: 2017-04-10 21:15:06 UTC
* Number of recursive dependencies: 73

Run `revdep_details(, "lspline")` for more info

</details>

## In both

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
      ...
    --- re-building ‘lspline.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    
    Quitting from lines 18-23 [setup] (lspline.Rmd)
    Error: processing vignette 'lspline.Rmd' failed with diagnostics:
    object 'params' not found
    --- failed re-building ‘lspline.Rmd’
    
    SUMMARY: processing the following file failed:
      ‘lspline.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# lucid

<details>

* Version: 1.8
* GitHub: https://github.com/kwstat/lucid
* Source code: https://github.com/cran/lucid
* Date/Publication: 2021-04-16 13:40:05 UTC
* Number of recursive dependencies: 73

Run `revdep_details(, "lucid")` for more info

</details>

## In both

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
      ...
    --- re-building ‘lucid_examples.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    
    Quitting from lines 234-258 [jags] (lucid_examples.Rmd)
    Error: processing vignette 'lucid_examples.Rmd' failed with diagnostics:
    could not find function "jags.model"
    --- failed re-building ‘lucid_examples.Rmd’
    
    SUMMARY: processing the following file failed:
      ‘lucid_examples.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘rjags’
    ```

# MagmaClustR

<details>

* Version: 1.1.2
* GitHub: https://github.com/ArthurLeroy/MagmaClustR
* Source code: https://github.com/cran/MagmaClustR
* Date/Publication: 2023-05-23 18:32:03 UTC
* Number of recursive dependencies: 108

Run `revdep_details(, "MagmaClustR")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘gifski’
    ```

# MapeBay

<details>

* Version: 0.1.0
* GitHub: NA
* Source code: https://github.com/cran/MapeBay
* Date/Publication: 2021-11-15 09:00:09 UTC
* Number of recursive dependencies: 188

Run `revdep_details(, "MapeBay")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 41.3Mb
      sub-directories of 1Mb or more:
        libs  41.0Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘RcppParallel’
      All declared Imports should be used.
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# mason

<details>

* Version: 0.3.0
* GitHub: https://github.com/lwjohnst86/mason
* Source code: https://github.com/cran/mason
* Date/Publication: 2020-06-04 16:10:05 UTC
* Number of recursive dependencies: 73

Run `revdep_details(, "mason")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# metabolomicsR

<details>

* Version: 1.0.0
* GitHub: https://github.com/XikunHan/metabolomicsR
* Source code: https://github.com/cran/metabolomicsR
* Date/Publication: 2022-04-29 07:40:02 UTC
* Number of recursive dependencies: 180

Run `revdep_details(, "metabolomicsR")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘genuMet’
    ```

# mhurdle

<details>

* Version: 1.3-0
* GitHub: NA
* Source code: https://github.com/cran/mhurdle
* Date/Publication: 2021-12-10 12:20:02 UTC
* Number of recursive dependencies: 76

Run `revdep_details(, "mhurdle")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

# mice

<details>

* Version: 3.15.0
* GitHub: https://github.com/amices/mice
* Source code: https://github.com/cran/mice
* Date/Publication: 2022-11-19 13:00:02 UTC
* Number of recursive dependencies: 134

Run `revdep_details(, "mice")` for more info

</details>

## Newly fixed

*   checking tests ...
    ```
      Running ‘testthat.R’
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      Backtrace:
          ▆
       1. ├─testthat::expect_warning(K <- parlmice(nhanes, n.core = 2, m = 7)) at test-parlmice.R:32:2
       2. │ └─testthat:::quasi_capture(...)
       3. │   ├─testthat (local) .capture(...)
       4. │   │ └─base::withCallingHandlers(...)
       5. │   └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
       6. └─mice::parlmice(nhanes, n.core = 2, m = 7)
       7.   └─parallel::makeCluster(n.core, type = cl.type)
       8.     └─parallel::makePSOCKcluster(names = spec, ...)
       9.       └─base::serverSocket(port = port)
      
      [ FAIL 5 | WARN 0 | SKIP 0 | PASS 367 ]
      Error: Test failures
      Execution halted
    ```

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

*   checking Rd \usage sections ... NOTE
    ```
    S3 methods shown with full name in documentation object 'cbind.mids':
      ‘cbind.mids’
    
    S3 methods shown with full name in documentation object 'rbind.mids':
      ‘rbind.mids’
    
    The \usage entries for S3 methods should use the \method markup and not
    their full name.
    See chapter ‘Writing R documentation files’ in the ‘Writing R
    Extensions’ manual.
    ```

# microbial

<details>

* Version: 0.0.20
* GitHub: NA
* Source code: https://github.com/cran/microbial
* Date/Publication: 2021-11-01 15:40:02 UTC
* Number of recursive dependencies: 172

Run `revdep_details(, "microbial")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘testthat’
      All declared Imports should be used.
    ```

# moderndive

<details>

* Version: 0.5.5
* GitHub: https://github.com/moderndive/moderndive
* Source code: https://github.com/cran/moderndive
* Date/Publication: 2022-12-01 22:20:02 UTC
* Number of recursive dependencies: 107

Run `revdep_details(, "moderndive")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 5238 marked UTF-8 strings
    ```

# MSEtool

<details>

* Version: 3.6.2
* GitHub: https://github.com/Blue-Matter/MSEtool
* Source code: https://github.com/cran/MSEtool
* Date/Publication: 2023-03-28 21:20:05 UTC
* Number of recursive dependencies: 147

Run `revdep_details(, "MSEtool")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 12.9Mb
      sub-directories of 1Mb or more:
        data   4.5Mb
        libs   5.6Mb
        R      1.9Mb
    ```

# multifear

<details>

* Version: 0.1.2
* GitHub: https://github.com/AngelosPsy/multifear
* Source code: https://github.com/cran/multifear
* Date/Publication: 2021-06-01 20:50:02 UTC
* Number of recursive dependencies: 119

Run `revdep_details(, "multifear")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘ez’
      All declared Imports should be used.
    ```

# multiverse

<details>

* Version: 0.6.1
* GitHub: https://github.com/MUCollective/multiverse
* Source code: https://github.com/cran/multiverse
* Date/Publication: 2022-07-04 13:20:02 UTC
* Number of recursive dependencies: 120

Run `revdep_details(, "multiverse")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘gifski’
    ```

# NetMix

<details>

* Version: 0.2.0.1
* GitHub: https://github.com/solivella/NetMix
* Source code: https://github.com/cran/NetMix
* Date/Publication: 2022-11-16 16:34:41 UTC
* Number of recursive dependencies: 59

Run `revdep_details(, "NetMix")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

# networkDynamic

<details>

* Version: 0.11.3
* GitHub: NA
* Source code: https://github.com/cran/networkDynamic
* Date/Publication: 2023-02-16 08:20:02 UTC
* Number of recursive dependencies: 39

Run `revdep_details(, "networkDynamic")` for more info

</details>

## In both

*   checking tests ...
    ```
      Running ‘activate_tests.R’
      Running ‘age.at_tests.R’
      Running ‘classTests.R’
      Running ‘converter_tests.R’
      Running ‘get_tests.R’
      Running ‘import_tests.R’
      Running ‘network_tests.R’
      Running ‘pid_tests.R’
      Running ‘query_tests.R’
      Running ‘reconcile.activityTests.R’
    ...
      > expect_equal(network.size(net),5)
      > expect_true(is.networkDynamic(net))
      > expect_equal(unlist(get.vertex.activity(net,as.spellList=TRUE)[4:5,1:2]),c(1,1,2,2),check.names=FALSE)
      Error: unlist(get.vertex.activity(net, as.spellList = TRUE)[4:5, 1:2]) (`actual`) not equal to c(1, 1, 2, 2) (`expected`).
      
      `names(actual)` is a character vector ('onset1', 'onset2', 'terminus1', 'terminus2')
      `names(expected)` is absent
      In addition: Warning message:
      Unused arguments (check.names = FALSE) 
      Execution halted
    ```

# nlmixr2est

<details>

* Version: 2.1.5
* GitHub: https://github.com/nlmixr2/nlmixr2est
* Source code: https://github.com/cran/nlmixr2est
* Date/Publication: 2023-04-22 19:50:02 UTC
* Number of recursive dependencies: 198

Run `revdep_details(, "nlmixr2est")` for more info

</details>

## In both

*   checking package dependencies ... ERROR
    ```
    Package required but not available: ‘symengine’
    
    See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
    manual.
    ```

# nlshelper

<details>

* Version: 0.2
* GitHub: NA
* Source code: https://github.com/cran/nlshelper
* Date/Publication: 2017-04-03 20:19:13 UTC
* Number of recursive dependencies: 38

Run `revdep_details(, "nlshelper")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# oceanic

<details>

* Version: 0.1.6
* GitHub: NA
* Source code: https://github.com/cran/oceanic
* Date/Publication: 2023-05-09 04:20:02 UTC
* Number of recursive dependencies: 53

Run `revdep_details(, "oceanic")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 1242 marked UTF-8 strings
    ```

# openintro

<details>

* Version: 2.4.0
* GitHub: https://github.com/OpenIntroStat/openintro
* Source code: https://github.com/cran/openintro
* Date/Publication: 2022-09-01 02:40:02 UTC
* Number of recursive dependencies: 93

Run `revdep_details(, "openintro")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  5.1Mb
      sub-directories of 1Mb or more:
        data   4.2Mb
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 2524 marked UTF-8 strings
    ```

# packDAMipd

<details>

* Version: 0.2.2
* GitHub: https://github.com/sheejamk/packDAMipd
* Source code: https://github.com/cran/packDAMipd
* Date/Publication: 2021-03-03 09:20:14 UTC
* Number of recursive dependencies: 217

Run `revdep_details(, "packDAMipd")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘flexsurv’ ‘nlme’ ‘tibble’ ‘tidyverse’ ‘tm’
      All declared Imports should be used.
    ```

# photosynthesis

<details>

* Version: 2.1.3
* GitHub: https://github.com/cdmuir/photosynthesis
* Source code: https://github.com/cran/photosynthesis
* Date/Publication: 2023-05-11 21:20:02 UTC
* Number of recursive dependencies: 143

Run `revdep_details(, "photosynthesis")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  7.2Mb
      sub-directories of 1Mb or more:
        doc   6.1Mb
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 13 marked UTF-8 strings
    ```

# pixiedust

<details>

* Version: 0.9.1
* GitHub: https://github.com/nutterb/pixiedust
* Source code: https://github.com/cran/pixiedust
* Date/Publication: 2021-01-15 11:50:02 UTC
* Number of recursive dependencies: 70

Run `revdep_details(, "pixiedust")` for more info

</details>

## In both

*   checking Rd files ... NOTE
    ```
    checkRd: (-1) sanitize_latex.Rd:49: Escaped LaTeX specials: \$
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# pkggraph

<details>

* Version: 0.2.3
* GitHub: https://github.com/talegari/pkggraph
* Source code: https://github.com/cran/pkggraph
* Date/Publication: 2018-11-15 09:50:03 UTC
* Number of recursive dependencies: 71

Run `revdep_details(, "pkggraph")` for more info

</details>

## In both

*   checking package subdirectories ... NOTE
    ```
    Problems with news in ‘NEWS.md’:
    No news entries found.
    ```

# plm

<details>

* Version: 2.6-3
* GitHub: https://github.com/ycroissant/plm
* Source code: https://github.com/cran/plm
* Date/Publication: 2023-04-09 11:40:02 UTC
* Number of recursive dependencies: 110

Run `revdep_details(, "plm")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  5.7Mb
      sub-directories of 1Mb or more:
        tests   3.1Mb
    ```

# plotly

<details>

* Version: 4.10.1
* GitHub: https://github.com/plotly/plotly.R
* Source code: https://github.com/cran/plotly
* Date/Publication: 2022-11-07 07:30:03 UTC
* Number of recursive dependencies: 167

Run `revdep_details(, "plotly")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  7.2Mb
      sub-directories of 1Mb or more:
        htmlwidgets   4.0Mb
    ```

# pubh

<details>

* Version: 1.2.7
* GitHub: https://github.com/josie-athens/pubh
* Source code: https://github.com/cran/pubh
* Date/Publication: 2022-04-04 13:50:02 UTC
* Number of recursive dependencies: 233

Run `revdep_details(, "pubh")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘Hmisc’ ‘sjPlot’
      All declared Imports should be used.
    ```

# qPCRtools

<details>

* Version: 0.2.1
* GitHub: https://github.com/lixiang117423/qPCRtools
* Source code: https://github.com/cran/qPCRtools
* Date/Publication: 2022-08-15 10:30:02 UTC
* Number of recursive dependencies: 137

Run `revdep_details(, "qPCRtools")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘data.table’ ‘readxl’ ‘reshape2’ ‘sjmisc’ ‘stringr’ ‘xlsx’
      All declared Imports should be used.
    ```

# radiant.data

<details>

* Version: 1.5.6
* GitHub: https://github.com/radiant-rstats/radiant.data
* Source code: https://github.com/cran/radiant.data
* Date/Publication: 2023-04-23 07:00:02 UTC
* Number of recursive dependencies: 146

Run `revdep_details(, "radiant.data")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘RPostgres’
    ```

# RBesT

<details>

* Version: 1.6-6
* GitHub: https://github.com/Novartis/RBesT
* Source code: https://github.com/cran/RBesT
* Date/Publication: 2023-03-03 18:20:02 UTC
* Number of recursive dependencies: 131

Run `revdep_details(, "RBesT")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 69.0Mb
      sub-directories of 1Mb or more:
        doc    1.8Mb
        libs  65.8Mb
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# rdss

<details>

* Version: 1.0.4
* GitHub: NA
* Source code: https://github.com/cran/rdss
* Date/Publication: 2023-05-02 08:10:03 UTC
* Number of recursive dependencies: 201

Run `revdep_details(, "rdss")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘DIDmultiplegt’
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 5 marked UTF-8 strings
    ```

# rfars

<details>

* Version: 0.3.0
* GitHub: https://github.com/s87jackson/rfars
* Source code: https://github.com/cran/rfars
* Date/Publication: 2023-05-05 09:40:02 UTC
* Number of recursive dependencies: 94

Run `revdep_details(, "rfars")` for more info

</details>

## In both

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 806 marked UTF-8 strings
    ```

# sand

<details>

* Version: 2.0.0
* GitHub: https://github.com/kolaczyk/sand
* Source code: https://github.com/cran/sand
* Date/Publication: 2020-07-02 07:20:06 UTC
* Number of recursive dependencies: 155

Run `revdep_details(, "sand")` for more info

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

# SherlockHolmes

<details>

* Version: 1.0.1
* GitHub: NA
* Source code: https://github.com/cran/SherlockHolmes
* Date/Publication: 2023-03-28 14:40:06 UTC
* Number of recursive dependencies: 108

Run `revdep_details(, "SherlockHolmes")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  6.4Mb
      sub-directories of 1Mb or more:
        doc       1.7Mb
        extdata   4.2Mb
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 14 marked UTF-8 strings
    ```

# simTool

<details>

* Version: 1.1.7
* GitHub: https://github.com/MarselScheer/simTool
* Source code: https://github.com/cran/simTool
* Date/Publication: 2020-09-22 16:00:03 UTC
* Number of recursive dependencies: 62

Run `revdep_details(, "simTool")` for more info

</details>

## Newly broken

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
      ...
    --- re-building ‘simTool.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    
    Quitting from lines 262-266 [unnamed-chunk-14] (simTool.Rmd)
    Error: processing vignette 'simTool.Rmd' failed with diagnostics:
    creation of server socket failed: port 11435 cannot be opened
    --- failed re-building ‘simTool.Rmd’
    
    SUMMARY: processing the following file failed:
      ‘simTool.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

# simts

<details>

* Version: 0.2.1
* GitHub: https://github.com/SMAC-Group/simts
* Source code: https://github.com/cran/simts
* Date/Publication: 2022-09-03 17:00:02 UTC
* Number of recursive dependencies: 55

Run `revdep_details(, "simts")` for more info

</details>

## In both

*   checking C++ specification ... NOTE
    ```
      Specified C++11: please drop specification unless essential
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 32.3Mb
      sub-directories of 1Mb or more:
        doc    1.5Mb
        libs  29.9Mb
    ```

# skedastic

<details>

* Version: 2.0.1
* GitHub: https://github.com/tjfarrar/skedastic
* Source code: https://github.com/cran/skedastic
* Date/Publication: 2022-11-06 07:40:02 UTC
* Number of recursive dependencies: 200

Run `revdep_details(, "skedastic")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘cmna’, ‘lmboot’
    ```

# sparklyr

<details>

* Version: 1.8.1
* GitHub: https://github.com/sparklyr/sparklyr
* Source code: https://github.com/cran/sparklyr
* Date/Publication: 2023-03-22 13:40:02 UTC
* Number of recursive dependencies: 116

Run `revdep_details(, "sparklyr")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  6.3Mb
      sub-directories of 1Mb or more:
        java   3.4Mb
        R      1.5Mb
    ```

# spqdep

<details>

* Version: 0.1.2
* GitHub: NA
* Source code: https://github.com/cran/spqdep
* Date/Publication: 2022-03-28 16:20:02 UTC
* Number of recursive dependencies: 108

Run `revdep_details(, "spqdep")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘lwgeom’ ‘rgeoda’
      All declared Imports should be used.
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 4 marked UTF-8 strings
    ```

# stabiliser

<details>

* Version: 1.0.6
* GitHub: NA
* Source code: https://github.com/cran/stabiliser
* Date/Publication: 2023-05-17 11:00:05 UTC
* Number of recursive dependencies: 149

Run `revdep_details(, "stabiliser")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘knitr’
      All declared Imports should be used.
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
      'brglm', 'dynlm', 'eha', 'erer', 'fGarch', 'gmm'
    ```

# statnet

<details>

* Version: 2019.6
* GitHub: https://github.com/statnet/statnet
* Source code: https://github.com/cran/statnet
* Date/Publication: 2019-06-14 08:00:06 UTC
* Number of recursive dependencies: 98

Run `revdep_details(, "statnet")` for more info

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

# statnet.common

<details>

* Version: 4.9.0
* GitHub: https://github.com/statnet/statnet.common
* Source code: https://github.com/cran/statnet.common
* Date/Publication: 2023-05-24 16:10:02 UTC
* Number of recursive dependencies: 20

Run `revdep_details(, "statnet.common")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    Found the following possibly unsafe calls:
    File ‘statnet.common/R/control.utilities.R’:
      unlockBinding("snctrl", environment(snctrl))
      unlockBinding("snctrl", environment(update_my_snctrl))
    ```

# statnetWeb

<details>

* Version: 0.5.6
* GitHub: NA
* Source code: https://github.com/cran/statnetWeb
* Date/Publication: 2020-08-05 18:00:03 UTC
* Number of recursive dependencies: 67

Run `revdep_details(, "statnetWeb")` for more info

</details>

## In both

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# StroupGLMM

<details>

* Version: 0.1.0
* GitHub: https://github.com/MYaseen208/StroupGLMM
* Source code: https://github.com/cran/StroupGLMM
* Date/Publication: 2016-04-19 01:00:56
* Number of recursive dependencies: 92

Run `revdep_details(, "StroupGLMM")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘broom’ ‘car’ ‘lmerTest’ ‘pbkrtest’
      All declared Imports should be used.
    ```

# survminer

<details>

* Version: 0.4.9
* GitHub: https://github.com/kassambara/survminer
* Source code: https://github.com/cran/survminer
* Date/Publication: 2021-03-09 09:50:03 UTC
* Number of recursive dependencies: 129

Run `revdep_details(, "survminer")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  6.0Mb
      sub-directories of 1Mb or more:
        doc   5.5Mb
    ```

# SWMPrExtension

<details>

* Version: 2.2.4.2
* GitHub: https://github.com/NOAA-OCM/SWMPrExtension
* Source code: https://github.com/cran/SWMPrExtension
* Date/Publication: 2023-04-20 22:12:34 UTC
* Number of recursive dependencies: 147

Run `revdep_details(, "SWMPrExtension")` for more info

</details>

## In both

*   checking whether package ‘SWMPrExtension’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/SWMPrExtension/new/SWMPrExtension.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘SWMPrExtension’ ...
** package ‘SWMPrExtension’ successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/srv/scratch/z3528859/R-4.3.0/library/magick/libs/magick.so':
  /apps/z_install_tree/linux-rocky8-ivybridge/gcc-12.2.0/glib-2.74.1-jrm3i3hfvjbyqoo4hipwrurprh3uteor/lib/libgobject-2.0.so.0: undefined symbol: g_uri_ref
Calls: <Anonymous> ... asNamespace -> loadNamespace -> library.dynam -> dyn.load
Execution halted
ERROR: lazy loading failed for package ‘SWMPrExtension’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/SWMPrExtension/new/SWMPrExtension.Rcheck/SWMPrExtension’


```
### CRAN

```
* installing *source* package ‘SWMPrExtension’ ...
** package ‘SWMPrExtension’ successfully unpacked and MD5 sums checked
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/srv/scratch/z3528859/R-4.3.0/library/magick/libs/magick.so':
  /apps/z_install_tree/linux-rocky8-ivybridge/gcc-12.2.0/glib-2.74.1-jrm3i3hfvjbyqoo4hipwrurprh3uteor/lib/libgobject-2.0.so.0: undefined symbol: g_uri_ref
Calls: <Anonymous> ... asNamespace -> loadNamespace -> library.dynam -> dyn.load
Execution halted
ERROR: lazy loading failed for package ‘SWMPrExtension’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/SWMPrExtension/old/SWMPrExtension.Rcheck/SWMPrExtension’


```
# SynDI

<details>

* Version: 0.1.0
* GitHub: https://github.com/umich-biostatistics/SynDI
* Source code: https://github.com/cran/SynDI
* Date/Publication: 2022-05-25 07:50:05 UTC
* Number of recursive dependencies: 66

Run `revdep_details(, "SynDI")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘arm’ ‘boot’ ‘broom’ ‘knitr’ ‘MASS’ ‘mvtnorm’ ‘randomForest’
      ‘StackImpute’
      All declared Imports should be used.
    ```

# takos

<details>

* Version: 0.2.0
* GitHub: https://github.com/sere3s/takos
* Source code: https://github.com/cran/takos
* Date/Publication: 2020-10-19 09:50:02 UTC
* Number of recursive dependencies: 45

Run `revdep_details(, "takos")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespaces in Imports field not imported from:
      ‘devEMF’ ‘MASS’ ‘segmented’ ‘smoother’ ‘tools’
      All declared Imports should be used.
    ```

*   checking LazyData ... NOTE
    ```
      'LazyData' is specified without a 'data' directory
    ```

# tergm

<details>

* Version: 4.1.1
* GitHub: https://github.com/statnet/tergm
* Source code: https://github.com/cran/tergm
* Date/Publication: 2022-11-08 14:10:02 UTC
* Number of recursive dependencies: 79

Run `revdep_details(, "tergm")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
      Running ‘degree.mean.age.R’
      Running ‘dynamic_EGMME.R’
      Running ‘dynamic_MLE_blockdiag.bipartite.R’
      Running ‘dynamic_MLE_blockdiag.R’
      Running ‘sim_gf_sum.R’
      Running ‘simulate_networkDynamic.R’
      Running ‘tergm_offset_tests.R’
      Running ‘tergm_parallel.R’
      Running ‘testthat.R’
     ERROR
    ...
        5.     └─base::replicate(...)
        6.       └─base::sapply(...)
        7.         └─base::lapply(X = X, FUN = FUN, ...)
        8.           └─tergm (local) FUN(X[[i]], ...)
        9.             └─tergm::tergm_MCMC_sample(...)
       10.               └─tergm::tergm_MCMC_slave(state, eta.comb, control, verbose)
      
      [ FAIL 1 | WARN 29 | SKIP 0 | PASS 3598 ]
      Error: Test failures
      Execution halted
    ```

# texPreview

<details>

* Version: 2.0.0
* GitHub: https://github.com/yonicd/texPreview
* Source code: https://github.com/cran/texPreview
* Date/Publication: 2022-03-31 07:30:02 UTC
* Number of recursive dependencies: 94

Run `revdep_details(, "texPreview")` for more info

</details>

## In both

*   checking tests ...
    ```
      Running ‘testthat.R’/srv/scratch/z3528859/R-4.3.0/bin/BATCH: line 60: 2818334 Aborted                 ${R_HOME}/bin/R -f ${in} ${opts} ${R_BATCH_OPTIONS} > ${out} 2>&1
    
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Complete output:
      > library(testthat)
      > library(texPreview)
      > test_check("texPreview")
      malloc(): invalid size (unsorted)
    ```

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
    sh: line 1: 2818815 Aborted                 '/srv/scratch/z3528859/R-4.3.0/bin/R' --vanilla --no-echo > '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d51238f27' 2>&1 < '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d53b769dde'
    --- re-building ‘classes.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    malloc(): invalid size (unsorted)
    sh: line 1: 2818881 Aborted                 '/srv/scratch/z3528859/R-4.3.0/bin/R' --vanilla --no-echo > '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d5467a47b' 2>&1 < '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d52a40d30e'
    --- re-building ‘engine.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    malloc(): invalid size (unsorted)
    sh: line 1: 2818939 Aborted                 '/srv/scratch/z3528859/R-4.3.0/bin/R' --vanilla --no-echo > '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d543786a0d' 2>&1 < '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d523a2decf'
    ...
    
    sh: line 1: 2819190 Aborted                 '/srv/scratch/z3528859/R-4.3.0/bin/R' --vanilla --no-echo > '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d539340f7' 2>&1 < '/scratch/pbs.4420393.kman.restech.unsw.edu.au/RtmpjxWJM3/file2b02d574e427e'
    --- re-building ‘tikz.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    malloc(): invalid size (unsorted)
    SUMMARY: processing the following files failed:
      ‘classes.Rmd’ ‘engine.Rmd’ ‘kable.Rmd’ ‘rmarkdown.Rmd’ ‘tikz.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘pdftools’
    ```

# texreg

<details>

* Version: 1.38.6
* GitHub: https://github.com/leifeld/texreg
* Source code: https://github.com/cran/texreg
* Date/Publication: 2022-04-06 22:00:02 UTC
* Number of recursive dependencies: 86

Run `revdep_details(, "texreg")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages which this enhances but not available for checking:
      'alpaca', 'brglm', 'dynlm', 'eha', 'erer', 'fGarch', 'gamlss.inf',
      'gmm', 'metaSEM', 'mnlogit', 'oglmx', 'pglm', 'simex', 'spatialreg',
      'Zelig'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘spatialreg’, ‘eha’, ‘brglm’, ‘dynlm’, ‘fGarch’, ‘alpaca’, ‘gamlss.inf’, ‘gmm’, ‘erer’, ‘oglmx’, ‘pglm’, ‘simex’, ‘metaSEM’
    ```

# tidybayes

<details>

* Version: 3.0.4
* GitHub: https://github.com/mjskay/tidybayes
* Source code: https://github.com/cran/tidybayes
* Date/Publication: 2023-03-14 04:30:02 UTC
* Number of recursive dependencies: 199

Run `revdep_details(, "tidybayes")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'runjags', 'rjags', 'jagsUI', 'gifski'
    ```

*   checking Rd cross-references ... NOTE
    ```
    Packages unavailable to check Rd xrefs: ‘runjags’, ‘rjags’, ‘jagsUI’
    ```

# tidycat

<details>

* Version: 0.1.2
* GitHub: https://github.com/guyabel/tidycat
* Source code: https://github.com/cran/tidycat
* Date/Publication: 2021-08-02 04:20:01 UTC
* Number of recursive dependencies: 71

Run `revdep_details(, "tidycat")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘tidyr’
      All declared Imports should be used.
    ```

# tidyquant

<details>

* Version: 1.0.7
* GitHub: https://github.com/business-science/tidyquant
* Source code: https://github.com/cran/tidyquant
* Date/Publication: 2023-03-31 20:40:02 UTC
* Number of recursive dependencies: 178

Run `revdep_details(, "tidyquant")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘tidyverse’
      All declared Imports should be used.
    ```

# timetk

<details>

* Version: 2.8.3
* GitHub: https://github.com/business-science/timetk
* Source code: https://github.com/cran/timetk
* Date/Publication: 2023-03-30 14:20:05 UTC
* Number of recursive dependencies: 177

Run `revdep_details(, "timetk")` for more info

</details>

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘timetk-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: parse_date2
    > ### Title: Fast, flexible date and datetime parsing
    > ### Aliases: parse_date2 parse_datetime2
    > 
    > ### ** Examples
    > 
    > 
    ...
    [1] "2011-01-01"
    > parse_date2("2011 June 3rd")
    [1] "2011-06-03"
    > 
    > # Fast datetime parsing
    > parse_datetime2("2011")
    Error in C_force_tz(to_posixct(time), tz, roll_dst) : 
      CCTZ: Unrecognized timezone of the input vector: ":/etc/localtime"
    Calls: parse_datetime2 ... <Anonymous> -> .force_tz -> from_posixct -> C_force_tz
    Execution halted
    ```

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
    --- re-building ‘TK04_Plotting_Time_Series.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    Warning in file.info(x, extra_cols = FALSE) :
      expanded path length 12200 would be too long for
    <p>This tutorial focuses on, <code>plot_time_series()</code>, a workhorse time-series plotting function that:</p>
    <ul>
    <li>Generates interactive <code>plotly</code> plots (great for exploring &amp; shiny apps)</li>
    <li>Consolidates 20+ lines of <code>ggplot2</code> &amp; <code>plotly</code> code</li>
    <li>Scales well to many time series</li>
    ...
    ℹ In group 1: `symbol = "AMZN"`.
    Caused by error in `C_time_add()`:
    ! CCTZ: Invalid timezone of the input vector: ""
    --- failed re-building ‘TK07_Time_Series_Data_Wrangling.Rmd’
    
    SUMMARY: processing the following file failed:
      ‘TK07_Time_Series_Data_Wrangling.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 2750 marked UTF-8 strings
    ```

# tipsae

<details>

* Version: 0.0.13
* GitHub: NA
* Source code: https://github.com/cran/tipsae
* Date/Publication: 2023-04-21 16:00:03 UTC
* Number of recursive dependencies: 156

Run `revdep_details(, "tipsae")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 66.3Mb
      sub-directories of 1Mb or more:
        doc    1.9Mb
        libs  62.3Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘sp’
      All declared Imports should be used.
    ```

*   checking data for non-ASCII characters ... NOTE
    ```
      Note: found 1 marked UTF-8 string
    ```

*   checking for GNU extensions in Makefiles ... NOTE
    ```
    GNU make is a SystemRequirements.
    ```

# TSS.RESTREND

<details>

* Version: 0.3.1
* GitHub: NA
* Source code: https://github.com/cran/TSS.RESTREND
* Date/Publication: 2020-08-02 20:10:02 UTC
* Number of recursive dependencies: 112

Run `revdep_details(, "TSS.RESTREND")` for more info

</details>

## In both

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘gimms’
    ```

# valr

<details>

* Version: 0.6.8
* GitHub: https://github.com/rnabioco/valr
* Source code: https://github.com/cran/valr
* Date/Publication: 2023-05-16 14:10:02 UTC
* Number of recursive dependencies: 173

Run `revdep_details(, "valr")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 14.2Mb
      sub-directories of 1Mb or more:
        libs  13.0Mb
    ```

# visR

<details>

* Version: 0.3.1
* GitHub: https://github.com/openpharma/visR
* Source code: https://github.com/cran/visR
* Date/Publication: 2022-08-17 22:10:03 UTC
* Number of recursive dependencies: 137

Run `revdep_details(, "visR")` for more info

</details>

## In both

*   checking re-building of vignette outputs ... ERROR
    ```
    Error(s) in re-building vignettes:
    --- re-building ‘CDISC_ADaM.Rmd’ using rmarkdown
    Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
    Warning in file.info(x, extra_cols = FALSE) :
      expanded path length 26884 would be too long for
    <h2 id="introduction">Introduction</h2>
    <p>This tutorial illustrates how a standard time-to-event analysis can be done very efficiently when the data set adheres to the <a href="https://www.cdisc.org/standards/foundational/adam/adam-basic-data-structure-bds-time-event-tte-analyses-v1-0">CDISC ADaM standard</a>. A more detailed time-to-event analysis with a more broad overview of visR’s functionality is presented in another vignette.</p>
    <pre><code class="language-r">library(ggplot2)
    library(visR)
    </code></pre>
    ...
    Quitting from lines 101-105 [table1_render_options_dt] (Time_to_event_analysis.Rmd)
    Error: processing vignette 'Time_to_event_analysis.Rmd' failed with diagnostics:
    invalid 'path' argument
    --- failed re-building ‘Time_to_event_analysis.Rmd’
    
    SUMMARY: processing the following file failed:
      ‘Time_to_event_analysis.Rmd’
    
    Error: Vignette re-building failed.
    Execution halted
    ```

# webSDM

<details>

* Version: 1.1-3
* GitHub: https://github.com/giopogg/webSDM
* Source code: https://github.com/cran/webSDM
* Date/Publication: 2023-03-14 13:50:02 UTC
* Number of recursive dependencies: 190

Run `revdep_details(, "webSDM")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘Matrix’
      All declared Imports should be used.
    ```

# whippr

<details>

* Version: 0.1.2
* GitHub: https://github.com/fmmattioni/whippr
* Source code: https://github.com/cran/whippr
* Date/Publication: 2022-09-09 07:00:02 UTC
* Number of recursive dependencies: 135

Run `revdep_details(, "whippr")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: ‘anomalize’
    ```

