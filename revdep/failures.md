# broom

<details>

* Version: 
* GitHub: https://github.com/statnet/ergm
* Source code: NA
* Number of recursive dependencies: 0

</details>

## Error before installation

### Devel

```






```
### CRAN

```






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
# ergmito

<details>

* Version: 0.3-0
* GitHub: https://github.com/muriteams/ergmito
* Source code: https://github.com/cran/ergmito
* Date/Publication: 2020-08-10 21:40:02 UTC
* Number of recursive dependencies: 63

Run `revdep_details(, "ergmito")` for more info

</details>

## In both

*   checking whether package ‘ergmito’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/ergmito/new/ergmito.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘ergmito’ ...
** package ‘ergmito’ successfully unpacked and MD5 sums checked
** using staged installation
checking whether the C++ compiler works... yes
checking for C++ compiler default output file name... a.out
checking for suffix of executables... 
checking whether we are cross compiling... no
checking for suffix of object files... o
checking whether we are using the GNU C++ compiler... yes
checking whether g++-10 accepts -g... yes
checking how to run the C++ preprocessor... g++-10 -E
checking for g++-10 option to support OpenMP... -fopenmp
configure: creating ./config.status
config.status: creating src/Makevars
** libs
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c RcppExports.cpp -o RcppExports.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c count_stats.cpp -o count_stats.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c ergmito_ptr.cpp -o ergmito_ptr.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c induced_submat.cpp -o induced_submat.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c network.cpp -o network.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c powersets.cpp -o powersets.o
g++ -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-z,relro -o ergmito.so RcppExports.o count_stats.o ergmito_ptr.o induced_submat.o network.o powersets.o -llapack -lblas -lgfortran -lm -lquadmath -fopenmp -L/usr/lib/R/lib -lR
installing to /home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/ergmito/new/ergmito.Rcheck/00LOCK-ergmito/00new/ergmito/libs
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) : 
  there is no package called ‘ergm’
Calls: <Anonymous> ... loadNamespace -> withRestarts -> withOneRestart -> doWithOneRestart
Execution halted
ERROR: lazy loading failed for package ‘ergmito’
* removing ‘/home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/ergmito/new/ergmito.Rcheck/ergmito’

```
### CRAN

```
* installing *source* package ‘ergmito’ ...
** package ‘ergmito’ successfully unpacked and MD5 sums checked
** using staged installation
checking whether the C++ compiler works... yes
checking for C++ compiler default output file name... a.out
checking for suffix of executables... 
checking whether we are cross compiling... no
checking for suffix of object files... o
checking whether we are using the GNU C++ compiler... yes
checking whether g++-10 accepts -g... yes
checking how to run the C++ preprocessor... g++-10 -E
checking for g++-10 option to support OpenMP... -fopenmp
configure: creating ./config.status
config.status: creating src/Makevars
** libs
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c RcppExports.cpp -o RcppExports.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c count_stats.cpp -o count_stats.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c ergmito_ptr.cpp -o ergmito_ptr.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c induced_submat.cpp -o induced_submat.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c network.cpp -o network.o
g++ -std=gnu++11 -I"/usr/share/R/include" -DNDEBUG  -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/Rcpp/include' -I'/home/pavel/Documents/Research/Software/statnet/ergm/revdep/library/ergmito/RcppArmadillo/include'   -fopenmp -DARMA_USE_OPENMP -DARMA_64BIT_WORD -fpic  -g -O2 -fdebug-prefix-map=/build/r-base-OT058M/r-base-4.0.2=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c powersets.cpp -o powersets.o
g++ -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-z,relro -o ergmito.so RcppExports.o count_stats.o ergmito_ptr.o induced_submat.o network.o powersets.o -llapack -lblas -lgfortran -lm -lquadmath -fopenmp -L/usr/lib/R/lib -lR
installing to /home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/ergmito/old/ergmito.Rcheck/00LOCK-ergmito/00new/ergmito/libs
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) : 
  there is no package called ‘ergm’
Calls: <Anonymous> ... loadNamespace -> withRestarts -> withOneRestart -> doWithOneRestart
Execution halted
ERROR: lazy loading failed for package ‘ergmito’
* removing ‘/home/pavel/Documents/Research/Software/statnet/ergm/revdep/checks/ergmito/old/ergmito.Rcheck/ergmito’

```
# fergm

<details>

* Version: 
* GitHub: https://github.com/statnet/ergm
* Source code: NA
* Number of recursive dependencies: 0

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
