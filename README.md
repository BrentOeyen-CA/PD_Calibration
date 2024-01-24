# On Probability of Default and its relation to Observed Default Frequency and a Common Factor

This repository contains some of the R codes used in the article [Oeyen, Celis Salazar , 2019](https://repository.uantwerpen.be/docman/irua/03b316/162930_2.pdf).

[Code](https://github.com/BrentOeyen-CA/PD_Calibration/tree/main/Codes/Eurostox_example.R) contains an extraction logic for EUROSTOXX data and how to use this data to create a realistic example of a ADF timeseries and how to extrapolate from this timeseries the asset correlation parameter and the Point-in-Time'ness parameter
![Picture Eurostoxx](https://github.com/BrentOeyen-CA/PD_Calibration/tree/main/Figures/Eurostox.pdf)

[Code](https://github.com/BrentOeyen-CA/PD_Calibration/tree/main/Codes/Simulation.R) contains the simulation logic for 5 scenarios and plots the impact on the estimation of the asset correlation parameter and the Point-in-Time'ness parameter using the proposed methodology of the auhors.
![ADF](https://github.com/BrentOeyen-CA/PD_Calibration/tree/main/Figures/ADF.pdf)
