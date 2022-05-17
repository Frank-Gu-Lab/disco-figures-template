# disco-figures-template

This repository is intended to be a template for generating custom figures using DISCO data.

## Environment Setup Steps:
1) Clone Repository to local machine
2) Create new python environment
3) Activate environment, navigate to this repository's home directory
4) Run `pip install requirements.txt` in terminal
5) Set up science plots matplotlib style
* For detailed instructions on installing the science plots matplotlib style and LaTex see: https://github.com/garrettj403/SciencePlots
6) (optional, for colour reproducibility) add disco library custom colour cycler file to your matplotlib styles

* i) find .matplotlib dir on your machine
`>>> import matplotlib as mpl`
`>>> mpl.get_configdir()`
* ii) show hidden files, navigate to dir/stylelib,
* iii) To use the default disco colours, move the discolib file associated with this repo:
    * `notebooks/utils/discolib.mplstyle` to 
    * `mpl.get_configdir()/stylelib`

* Colours are from the colorbrewer qualitative comparison palette: https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9