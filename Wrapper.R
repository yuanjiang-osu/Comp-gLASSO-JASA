### Simulations for Dense Data (Figures 1-2 in Section 3.2)

setwd("./Simulation")
source("Simulation_Dense.R")

### Simulations for Sparse Data (Figures S1-S2 in Supplementary Material)

setwd("./Simulation")
source("Simulation_Sparse.R")

### TARA Data Analysis (Section 4.1)

setwd("./TARA")

# Generate intermediate results
# Comment to save time by using existing intermediate results
# source("data_preprocessing.R")
# source("TARA_path.R")
# source("TARA_StARS.R")

# Figure 3(a)
source("TARA_table_plot_path.R")

# Figures S3 and 3(b) and Tables 1 and S1
source("TARA_table_plot_StARS.R")

### Zebrafish Data Analysis (Section 4.1)

setwd("./Zebrafish")

# Generate intermediate results
# Comment to save time by using existing intermediate results
# source("data_preprocessing.R")
# source("fish_StARS.R")

# Figures 4-5 and Table S2
source("fish_table_plot_StARS.R")

