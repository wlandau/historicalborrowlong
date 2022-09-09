
# historicalborrowlong

[![check](https://github.com/wlandau/historicalborrowlong/workflows/check/badge.svg)](https://github.com/wlandau/historicalborrowlong/actions?query=workflow%3Acheck)
[![cover](https://github.com/wlandau/historicalborrowlong/workflows/cover/badge.svg)](https://github.com/wlandau/historicalborrowlong/actions?query=workflow%3Acover)
[![pkgdown](https://github.com/wlandau/historicalborrowlong/workflows/pkgdown/badge.svg)](https://github.com/wlandau/historicalborrowlong/actions?query=workflow%3Apkgdown)
[![lint](https://github.com/wlandau/historicalborrowlong/workflows/pkgdown/badge.svg)](https://github.com/wlandau/historicalborrowlong/actions?query=workflow%3Alint)

Historical borrowing in clinical trials can improve precision and
operating characteristics. This package supports a longitudinal
hierarchical model to borrow historical control data from other studies
to better characterize the control response of the current study. It
also quantifies the amount of borrowing through longitudinal benchmark
models (independent and pooled). The hierarchical model approach to
historical borrowing is discussed in [Viele et
al. (2013)](https://onlinelibrary.wiley.com/doi/10.1002/pst.1589).

## Installation

``` r
remotes::install_github("wlandau/historicalborrowlong")
```

## Documentation

-   Functions:
    <https://wlandau.github.io/historicalborrowlong/reference/>
-   Methods:
    <https://wlandau.github.io/historicalborrowlong/articles/methods.html>
-   Usage:
    <https://wlandau.github.io/historicalborrowlong/articles/usage.html>

## Thanks

[Albert Man](https://github.com/albert-man), [Faith
Bian](https://github.com/faithbian-lilly), and [Saptarshi
Chatterjee](https://github.com/schatterjee-lilly) contributed to the
development of the methods prior to implementation. Phebe Kemmer,
Heather Zhao, and Zhangchen Zhao also provided helpful feedback on the
models and their application to clinical use-cases.

## Code of Conduct

Please note that the historicalborrowlong project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.