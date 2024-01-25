# On Probability of Default and its relation to Observed Default Frequency and a Common Factor

This repository contains some of the R codes used in the article Oeyen, Brent & Salazar Celis, Oliver, **On probability of default and its relation to observed default frequency and a common factor**, *Journal of Credit Risk*, 15(3):41-66, 2019 ([doi](https://doi.org/10.21314/JCR.2019.253)).

[Code](https://github.com/BrentOeyen-CA/PD_Calibration/tree/main/Codes/Eurostox_example.R) contains an extraction logic for EUROSTOXX data and how to use this data to create a realistic example of a ADF timeseries and how to extrapolate from this timeseries the asset correlation parameter and the Point-in-Time'ness parameter
![eurostocks](https://github.com/BrentOeyen-CA/PD_Calibration/assets/7952417/064be11f-7b1b-42a6-b417-751656380260)

[Code](https://github.com/BrentOeyen-CA/PD_Calibration/tree/main/Codes/Simulation.R) contains the simulation logic for 5 scenarios and plots the impact on the estimation of the asset correlation parameter and the Point-in-Time'ness parameter using the proposed methodology of the auhors.
![adf](https://github.com/BrentOeyen-CA/PD_Calibration/assets/7952417/ceadd12c-f5c9-4565-8dd0-cdba1c4426b7)
