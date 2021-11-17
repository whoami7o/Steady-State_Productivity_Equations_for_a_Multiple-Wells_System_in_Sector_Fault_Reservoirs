# Steady-State Productivity Equations for a Multiple-Wells System in Sector Fault Reservoirs
This module created for use of findings presented in the following SPE article published by *"Jing Lu, Shawket Ghedan, and Tao Zhu, The Petroleum Institute; and Anh Dinh, Schlumberger; and Djebbar Tiab, University of Oklahoma"*.
For questions regarding to equations, their implementation or calc. procedure or more detailed information please refer to the article (it can be found in repository).
## Function library 
- **ssp_eq_mws_sfr.py**
Python file which implements all equations from article as functions for further use. They were tested on data and results presented in the article.
There are possibility for implementing field-units, but it can be done in the future.
## Calculation app
- **calculaion_app_gui.py**
Deploy this file to do your calculation. Interface is made as much as user-friendly as possible, so I hope you won't have any struggles working with it.
By default all cells are filled with values for **"Example One, *Case 1*"** as in the SPE article.
### Buttons
- **find Q** -  enables calculation mode when *skin-factor* of each well is known and outputs production rate ($Q, m^3/d$) for i-*th* well and their sum rate in the corresponding cell;
- **find S** -  enables calculation mode when production rate ($Q, m^3/d$) of each well is known and outputs *skin-factor* for i-*th* well in the corresponding cell;
- **+, -** - add/remove well;
- **Calculate** -  does what it says and outputs result;
- **Save results** - creates Excel *'manual_save_test.xlsx'* file with all your entered and calculated parameters in the same directory where script was deployed;
- **Calculate from excel** - runs *'cmd_only_excel_mode.py'* functionality and saves results in *'output_solution_example.xlsx'* in the same directory where script was deployed.
### Errors
All common errors  were catched preventing app from crashing. So, when one of them occurs the following self-explanatory message prints out in the console.
### Variables explanation
- $P_e$ (P_e) – reservoir outer-boundary pressure, $MPa$;
- $Φ$ (PHI) – angle of the sector reservoir, *boundaries* : $[0, 360]$, $° (deg)$;
- $Q$ (Q) - production rate, $m^3/d$;
- $r_e$ (r_e) - radius of the sector-fault reservoir, $m$;
- $H$ - pay-zone thickness, $m$;
- $K_r$ (K_r) – radial permeability,  $μm^2$;
- $K_z$ (K_z) – vertical permeability, $μm^2$;
- $μ$ (mu) – oil viscosity, $mPa·s$;
- $B$ – formation volume factor, $m^3/m^3$;
- $p_wf$ (P_wf) – well-flow pressure, $MPa$;
- $r_i$ (r_i) – off-vertex distance of i-th well in multiple-wells system, $m$;
- $φ$ (phi_i) – wellbore location angle of jth well in the multiple wells system, $° (deg)$;
- $r_w$ (r_w) – wellbore radius, $m$;
- $S$ (skin) – skin-factor of the well.

## CMD mode for Excel
- **cmd_only_excel_mode.py**
For those who don't want to deal with GUI application, calculation from Excel can be done by simply running this script. Be careful with it, especially with file path (I personally recommend to specify only absolute paths to input excel files).
The results will be saved into *'cmd_solution.xlsx'* in the same directory where script has been deployed.
## Input/output templates for Excel
### Input
- **input_data_template.xlsx**
The excel input allows you to run as much cases as you want without entering it manually. But there are some **rules** to creating input excel file (you can find template in the following file).
It's important to keep structure **exactly** as it is in template, except you can set as many well as you need and you can specify your input case by changing 'SKIN' to 'Q' depending on your known data.
(!CHANGE ONLY VALUES and NO BLANK ROWS ALLOWED IN CASE SCENARIO (between 'START and 'END'!))
### Output
File created by different 'saving mode' (be careful because their name doesn't depend on how many time you have run the scripts results will be saved in file with same name each time, so it can be re-written):
- **output_solution_example.xlsx**
- **manual_save_test.xlsx**
- **cmd_solution.xlsx**

## Extras
- **images** folder
Images for GUI app;
- **calculation_app_gui.ui**
PyQt5 Designer file.
- **requirements.txt**
