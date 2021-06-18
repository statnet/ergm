# btergm

<details>

* Version: 1.9.13
* GitHub: https://github.com/leifeld/btergm
* Source code: https://github.com/cran/btergm
* Date/Publication: 2020-10-26 14:30:02 UTC
* Number of recursive dependencies: 75

Run `revdep_details(, "btergm")` for more info

</details>

## In both

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
Error: object ‘ergm.Cprepare’ is not exported by 'namespace:ergm'
Execution halted
ERROR: lazy loading failed for package ‘btergm’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/btergm/old/btergm.Rcheck/btergm’


```
# tergm

<details>

* Version: 3.7.0
* GitHub: https://github.com/statnet/tergm
* Source code: https://github.com/cran/tergm
* Date/Publication: 2020-10-15 12:40:03 UTC
* Number of recursive dependencies: 46

Run `revdep_details(, "tergm")` for more info

</details>

## In both

*   checking whether package ‘tergm’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/tergm/new/tergm.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘tergm’ ...
** package ‘tergm’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
gcc -I"/srv/scratch/z3528859/R-4.1.0/include" -DNDEBUG  -I'/srv/scratch/z3528859/R-4.1.0/library/ergm/include' -I/usr/local/include   -fpic  -g -O2  -c DynSA.c -o DynSA.o
In file included from MCMCDyn.h:13,
                 from DynSA.h:13,
                 from DynSA.c:10:
/srv/scratch/z3528859/R-4.1.0/library/ergm/include/edgetree.h:10:9: note: #pragma message: warning: The header file "edgetree.h" has been deprecated in favor of "ergm_edgetree.h" and may be removed in the future.
 #pragma message ("warning: The header file \"edgetree.h\" has been deprecated in favor of \"ergm_edgetree.h\" and may be removed in the future.")
...
                ^~
DynSA.c:92:11: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     if(nwp->duration_info.lasttoggle)
           ^~
DynSA.c:93:27: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     memcpy(lasttoggle, nwp->duration_info.lasttoggle, sizeof(int)*DYADCOUNT(*n_nodes, *bipartite, *dflag));
                           ^~
make: *** [DynSA.o] Error 1
ERROR: compilation failed for package ‘tergm’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/tergm/new/tergm.Rcheck/tergm’


```
### CRAN

```
* installing *source* package ‘tergm’ ...
** package ‘tergm’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
gcc -I"/srv/scratch/z3528859/R-4.1.0/include" -DNDEBUG  -I'/srv/scratch/z3528859/R-4.1.0/library/ergm/include' -I/usr/local/include   -fpic  -g -O2  -c DynSA.c -o DynSA.o
In file included from MCMCDyn.h:13,
                 from DynSA.h:13,
                 from DynSA.c:10:
/srv/scratch/z3528859/R-4.1.0/library/ergm/include/edgetree.h:10:9: note: #pragma message: warning: The header file "edgetree.h" has been deprecated in favor of "ergm_edgetree.h" and may be removed in the future.
 #pragma message ("warning: The header file \"edgetree.h\" has been deprecated in favor of \"ergm_edgetree.h\" and may be removed in the future.")
...
                ^~
DynSA.c:92:11: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     if(nwp->duration_info.lasttoggle)
           ^~
DynSA.c:93:27: error: ‘Network’ {aka ‘struct Networkstruct’} has no member named ‘duration_info’
     memcpy(lasttoggle, nwp->duration_info.lasttoggle, sizeof(int)*DYADCOUNT(*n_nodes, *bipartite, *dflag));
                           ^~
make: *** [DynSA.o] Error 1
ERROR: compilation failed for package ‘tergm’
* removing ‘/srv/scratch/z3528859/github/statnet/ergm/revdep/checks/tergm/old/tergm.Rcheck/tergm’


```
