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

# btergm

<details>

* Version: 1.9.13
* GitHub: https://github.com/leifeld/btergm
* Source code: https://github.com/cran/btergm
* Date/Publication: 2020-10-26 14:30:02 UTC
* Number of recursive dependencies: 71

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
Warning message:
replacing previous import ‘vctrs::data_frame’ by ‘tibble::data_frame’ when loading ‘dplyr’ 
Error: object ‘ergm.Cprepare’ is not exported by 'namespace:ergm'
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
Warning message:
replacing previous import ‘vctrs::data_frame’ by ‘tibble::data_frame’ when loading ‘dplyr’ 
Warning: no DISPLAY variable so Tk is not available
** help
...
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
Warning: replacing previous import ‘vctrs::data_frame’ by ‘tibble::data_frame’ when loading ‘dplyr’
Warning: no DISPLAY variable so Tk is not available
** testing if installed package can be loaded from final location
Warning: replacing previous import ‘vctrs::data_frame’ by ‘tibble::data_frame’ when loading ‘dplyr’
Warning: no DISPLAY variable so Tk is not available
** testing if installed package keeps a record of temporary installation path
* DONE (btergm)


```
# EpiModel

<details>

* Version: 2.0.3
* GitHub: https://github.com/statnet/EpiModel
* Source code: https://github.com/cran/EpiModel
* Date/Publication: 2020-11-09 21:40:13 UTC
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
    Finished MPLE.
    Stopping at the initial estimate.
    > 
    > # Simulate the network with netdx
    > dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
    +             verbose = FALSE)
    Error in get(name, pos = environment(new)) : 
      argument "set.control.stergm" is missing, with no default
    Calls: netdx ... .handle.auto.constraints -> nonsimp_update.formula -> assign -> get
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/test-all.R’ failed.
    Last 13 lines of output:
       13.             │   └─base::eval(cl, parent.frame())
       14.             ├─stats::simulate(...)
       15.             └─ergm:::simulate.formula_lhs_network(...)
       16.               ├─ergm::simulate_formula(...)
       17.               └─ergm:::simulate_formula.network(...)
       18.                 └─ergm:::.handle.auto.constraints(...)
       19.                   └─statnet.common::nonsimp_update.formula(...)
       20.                     ├─base::assign(name, get(name, pos = environment(new)), pos = e)
       21.                     └─base::get(name, pos = environment(new))
      
      [ FAIL 33 | WARN 0 | SKIP 78 | PASS 231 ]
      Error: Test failures
      In addition: Warning message:
      replacing previous import 'vctrs::data_frame' by 'tibble::data_frame' when loading 'dplyr' 
      Execution halted
    ```

