# HZg_scripts
Various scripts used for CMS Higgs to Z gamma analysis

To apply BDT score to pico samples, use producer_Input.cxx. You will need to change the xml file location for the BDT and input variables. A lot of the jet variables are not currently saved to the reduced pico with BDT the producer makes. See below for an example command. You can make a .sh script to run all the years and samples in parallel in batch. 
```
root -b -q 'producer_Input.cxx("GGF", "Name of the output dir (need to be created beforehand)", "Name of the input dir with BDT score", "2016")'
```

Once you have an ntuple with BDT score saved in a branch. You can run optimization tests on it and make plots. 
To run an optimization test, use the following command:
```
root -b -q 'optbins.cxx(false (plot option, false means no plotting of the bdt score and mllg of the best found binning), true (true means use test sample only), ""(name tag of the output file, only used when plot option is true))'
```

To run plotting scripts with certain BDT score binning strategies, use the following command:
```
root -b -q 'plot_bdt.cxx(true, "", true, true, true, true)'
```
More plotting scripts can be uploaded.  

Both optbins.cxx and plot_bdt.cxx use RDataFrame, which is only introduced after ROOT version 6.14. For UCSB cluster, it means one needs to be on cms36.  
