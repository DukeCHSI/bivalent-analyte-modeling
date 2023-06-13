## Bivalent analyte package

This provides information about R-package for model fitting binding kinetics data using bivalent analyte.

## Folder Structure

The project comprises of two main folders:

### data/

This folder should contain two important files:

1. **<...> - info sheet.csv**: this csv file is an info sheet file that has the information about the data sets.
  - The package can handle these options: 
    - Regeneration: 
      - `Regen. = N` (default): NO regenerative titration.
      - `Regen. = Y`: regenerative titration.
    - Global or local Rmax:
      - `Global Rmax = N` (default): local Rmax, fit Rmax separately.
      - `Global Rmax = Y`: global Rmax, only fit one Rmax for each sensorgram.
    - Second rebinding:
      - `Sec. Rebinding = N` (default): assume that the second arms cannot rebind during dissociation phase.
      - `Sec. Rebinding = Y`: assume that the second arms can rebind during dissociation phase.
  - Other notes:
    - `Bulkshift`: is (possibly) available with another 1:1 package (more user-friendly), but currently not implemented in this package.
    - `Model`: this option was included because there is a separated version that could handle different models other than just the 1:1 model and the bivalent analyte model. However, it was developed on parallel with the more user-friendly version of the 1:1 model package. So, the implementation to handle other models have been removed.
    - `Refractive index`: the current version can handle refractive index. However, it is not an available option in the info sheet. If the user want to use refractive index, they could add another column called `Refractive index` to the info sheet file (.csv). Then, go to **ode_function.R** and change line 128 to:
      - use_RI <- sample_info[well_idx,]$`Refractive index`
2. **<...> - new output format.xlsx**: this xlsx file is an output file has the Time and Response Unit data values.

### code/

This folder should contain `.R` files and subfolders

In this folder:

- **helper_functions.R** contains most of the non-ODE related functions for data processing, data plotting.
- **ode_functions.R** contains ODE-related functions for solving the bivalent analyte model and plotting the results.
- **main.R** is the main script to run bivalent analyte fitting.
- **fitting_results/** is the folder that contains fitting result files (`.Rdata`).
- **figures/** contains the figure files

## Guide 

### Start guide

To start the bivalent analyte model fitting:

1. Make sure your info sheet and output files have the same format as the provided example files. 
2. Run the script **main.R**.
3. The script will ask to set the working directory. Choose the **code/** as your working directory.
4. The script will ask for the info sheet file.
5. The script will ask for the output file.

### Results

After fitting running,

1. The fitting results are saved in an `.Rdata` file within **fitting_results/**.
2. The fitting figures and sensorgrams are saved in **figures/bivalent/**.

### Troubleshooting

The required formatting for this package is not the same as other packages. If you run into problems,

1. In the info sheet file (.csv), make sure the numbers in `Block/Chip/Tray` are separated by semicolon (`;`).
2. In the output file (.xlsx), make sure there are 3 header rows. One possible option is to copy and inserted the header rows from the provided output file to your own output file.
3. In the info sheet file (.csv), make sure you have two columns `Sec. Rebinding` and `Model`. 