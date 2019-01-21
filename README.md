# Diverse Stochasticity Leads a Colony of Ants to Optimal Foraging

A foraging ants multi-agent simulation software used in the following publication.

## Publication
- Masashi Shiraishi, Rito Takeuchi, Hiroyuki Nakagawa, Shin I Nishimura, Akinori Awazu, and Hiraku Nishimori, "Diverse Stochasticity Leads a Colony of Ants to Optimal Foraging", [Journal of Theoretical Biology](https://doi.org/10.1016/j.jtbi.2019.01.002), 2019.

## Requirements
- Ubuntu >= 17.04 or Windows >= 7
- CUDA 8.0 or higher
    
    This software is initially designed under CUDA 5.5 and maybe works with this.

- freeglut
- Python 3
- GNU make

## Usage
```bash
[Ant-Cuda] python3 run.py -h
usage: run.py [-h] --angle ANGLE --dist DIST

A foraging ants multi-agent simulation software. The result is output into
subdirectories.

optional arguments:
  -h, --help     show this help message and exit
  --angle ANGLE  Relative angle between two food resources $\theta$.
  --dist DIST    Distance between the nest and each food $R$.
```

## Parameters
Parameters are defined at `Constants.h`.

* `MACRO_NMAX`: # of ants ![](https://latex.codecogs.com/gif.latex?N_{total}).
* `MACRO_MAX`: # of horizontal/vertical cells.
* `MACRO_FOODSOURCE`: Initial amount of each food resource.
* `MACRO_NUM_FOODS`: # of food resource sites.
* `MACRO_FOOD_DIST`: Distance between the nest and each food ![](https://latex.codecogs.com/gif.latex?R).
* `MACRO_FOOD_ANGLE`: Relative angle between two food resources ![](https://latex.codecogs.com/gif.latex?\theta).
* `MACRO_NEST_X`: Vertical position of the nest cell.
* `MACRO_NEST_Y`: Horizontal position of the nest cell.
* `MACRO_REC`: Recovery Ratio of each food resource ![](https://latex.codecogs.com/gif.latex?r).
* `MACRO_EVAPOLATION_CONST`: Evaporation rate of pheromone ![](https://latex.codecogs.com/gif.latex?e_p).
* `MACRO_MAX_SEARCH_TIME`: Time interval after which an ant is in Exploring mode ![](https://latex.codecogs.com/gif.latex?t_{th}).
* `MACRO_MAX_STEP`: # of ensemble of Monte-Carlo simulation for each ![](https://latex.codecogs.com/gif.latex?N_s) and ![](https://latex.codecogs.com/gif.latex?\alpha_s).
* `MACRO_MAX_TIME`: Max time step at each trial ![](https://latex.codecogs.com/gif.latex?T).
* `MACRO_EMI`: Secreted pheromone amount ![](https://latex.codecogs.com/gif.latex?\phi_u).
* `MACRO_UNIT`: Diminution of food resource at each foraging.
* `MACRO_DIFFE`: Diffusion constant ![](https://latex.codecogs.com/gif.latex?D).
