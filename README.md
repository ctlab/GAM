# GAM: integrated analysis of transcriptional and metabolic profiling

GAM consists of three packages: 
* Main `GAM` package that contains algorithms for network creation and integrated analysis
* `GAM.db` package that contains processed data from KEGG and HMDB
* `GAM.networks` package that containts networks for human and mouse.

First, run the `install-dependencies.R` script to install some necessary packages:
```
Rscript install-dependencies.R
```

Please, check that everything was install correctly during the above command!

Now you need to build and install `GAM.db` package:
* Run `make GAM.db-package` on the machine with access to KEGG to build the package.
* Run `R CMD INSTALL $(ls GAM.db_*.tar.gz | tail -1)` to install the package.

Second, you need to build `GAM.networks` package:
* Run `make GAM.networks-package`.
* Run `R CMD INSTALL $(ls GAM.networks_*.tar.gz | tail -1)` to install the package.

Finally, you can build the main `GAM` package by calling `make GAM`.

You can consult with `GAM/vignettes/GAM-tutorial.Rmd` on how to run the analysis.

## Solvers

GAM requires an MWCS solver to be available. For Linux x64 machines we recommend to download heinz (http://homepages.cwi.nl/~klau/data/heinz_1.68.tgz) and unpack it into `/usr/local/lib/heinz` directory.
