
# The Energy Cascade(s)

This project contains implementations of the Cavazzoni 2004, Boscheri 2012, Amitrano 2020, and Volk 1995 crop growth models. Previously theses were only available in outdated properitary formats or only as equations in their respective publications. If you are interested in a nitrogen productivity version of the Cavazzoni model, visit [kmyates262](https://github.com/kmyates262/nitrogen-productivity).

These models work of the approach of cascading photonic energy through plant physiological processes to predict crop yield. Originally developed for advanced life support systems these models have been developed since 1995 and are highly suitable for controlled environment agriculture.

## Overview of the Energy Cascade(s)

At its highest levels, the Energy Cascade(s) model photon absorption, canopy quantum yield, and carbon use efficiency. With that information they then calculate plant processes such as crop growth rate, photosynthesis rates, and gas exchanges to predict outputs of interest such as crop biomass, CO2 consumption, O2 production, and transpiration rate.

## Parameterization Note

The Energy Cascade(s) exist in various degrees of parameterization.

* Volk et. al was only parameterized for wheat
* Cavazzoni et. al was parameterized for soybeans
  * Provided untested parameters for dry bean, lettuce, peanut, rice, sweet potato, tomato, wheat, white potato
* Boscheri et. al used the Cavazzoni parameters
* Amitrano et. al was paramaterized for red and green lettuce

This project is currently only designed to work with lettuce. Future iterations will include the other parameters provided by Cavazzoni et. al.

## Features

* Implements Cavazzoni et al. 2004, Boscheri et. al 2012, Amitrano et al. 2020, and Volk et. al 1995 crop growth models.
* Reads daily environment and crop data from CSV files.
* Robust input validation and flexible output management.
* Cross-platform, open-source, and reproducible for research applications.

## Getting Started

### Prerequisites

* Conda (recommended for environment management)
* Python 3.10

### Installation

Clone the repository and set up the environment using the provided environment.yml. If you do not have Conda installed, please refer to the Conda installation guide:

```
git clone https://github.com/yourusername/crop-simulation.git
cd crop-simulation
conda env create -f environment.yml
conda activate crop-simulation
```

This will create and activate a new environment with all necessary packages for running the crop simulation code. Make sure that any local modules (Init.py, Dynamic_ECs.py) are present in your working directory.

## Usage

Review and adapt sample demonstration scripts (see main.py or similar) or import the models directly in your code:

```
from Dynamic_ECs import CAV, BOS, AMI, EC

results = CAV(input_file="input_file.csv", TCB=10, t_A=28, output="results", save_output=True)
```
See demonstration.py for an example of how to run the models

### Input Data Requirements

Input CSV files must contain columns the columns: t, TEMP, RH, PPFD, CO2, H, and P_ATM where

* t = timestep (days)
* TEMP = Daily Average Temperature (Celsius)
* RH = Daily Average Relative Humidity (%)
* PPFD = Daily Average Photosyntetic Photon Flux Density (µmol m^-2 sec^-1)
* CO2 = Daily Average Carbon Dioxide Concentration (ppm)
* H = Photoperiod (Hours)
* P_ATM = Daily Average Atmospheric Pressure (kPa)

See provided examples in the `input_file.csv` or consult the function docstrings for formatting.

### Output

Results are saved as CSV files in the designated results directory if save_output=True.

## Project Structure

* Init.py – Initialization routines, parameter management
* Dynamic_ECs.py – Core crop model implementations (CAV, BOS, AMI, EC)
* environment.yml – Conda environment file
* README.md – Project documentation
* input_file.csv – Example input data (please provide or reference)

## Contributing

Contributions and improvements are welcome—please fork the repository and submit pull requests!

## License

Distributed under the MIT License. See LICENSE for details.

## References

1. Cavazzoni, J. Using Explanatory Crop Models to Develop Simple Tools for Advanced Life Support System Studies. *Advances in Space Research*, 2004, 34, 1528–1538. [https://doi.org/10.1016/j.asr.2003.02.073](https://doi.org/10.1016/j.asr.2003.02.073)

2. Boscheri, G.; Kacira, M.; Patterson, L.; Giacomelli, G.; Sadler, P.; Furfaro, R.; Lobascio, C.; Lamantea, M.; Grizzaffi, L. Modified Energy Cascade Model Adapted for a Multicrop Lunar Greenhouse Prototype. *Advances in Space Research*, 2012, 50, 941–951. [https://doi.org/10.1016/j.asr.2012.05.025](https://doi.org/10.1016/j.asr.2012.05.025)

3. Amitrano, C.; Chirico, G.B.; De Pascale, S.; Rouphael, Y.; De Micco, V. Crop Management in Controlled Environment Agriculture (CEA) Systems Using Predictive Mathematical Models. *Sensors*, 2020, 20, 3110. [https://doi.org/10.3390/s20113110](https://doi.org/10.3390/s20113110)

4. Volk, T.; Bugbee, B.; Wheeler, R.M. An Approach to Crop Modeling with the Energy Cascade. *Life Support Biosph Sci*, 1995, 1, 119–127. [PMID: 11538584](https://pubmed.ncbi.nlm.nih.gov/11538584)

5. M. K. Ewert, T. T. Chen, and C. D. Powell, Life Support Baseline Values and Assumptions Document, Life Support, p. 235. [NASA Technical Report](https://ntrs.nasa.gov/api/citations/20180001338/downloads/20180001338.pdf)
