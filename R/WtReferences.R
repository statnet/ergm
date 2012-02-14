# For now, this file contains information about reference
# meausres. Eventually, we should create an "InitReference" or similar
# framework.


# A lookup table for methods for initial fits; the default is the first one.
init.methods<-list(
                   Bernoulli=c("MPLE","zeros"),
                   Poisson=c("zeros"),
                   DescRank=c("zeros"),
                   StdNormal=c("zeros")
                   )
