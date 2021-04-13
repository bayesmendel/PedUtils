# PedUtils
simulations+model evaluations

Dependencies: `PanelPRO`, `abind`

Notes from Jane 4/13/2021
- I added my family simulation functions this repo (prefixed with `sim.`). They are compatible with (and dependent on) the current PanelPRO package. 
- Since the PanelPRO package is still being updated, the simulation functions may also need updates. Please let me know if you run into any issues. 
- Example usage is below; it calls the sim.runSimFam wrapper function for the entire simulation process: 

```
# Cancers
cancers = c("Brain", "Breast", "Colorectal", "Endometrial", 
            "Gastric", "Kidney", "Melanoma", "Ovarian", 
            "Pancreas", "Prostate", "Small Intestine")
# Genes
genes = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM", 
          "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")
          
# Paternal aunts, paternal uncles
nSibsPatern = c(1, 2) 
# Maternal aunts, maternal uncles
nSibsMatern = c(2, 1) 
# Sisters and brothers
nSibs = c(2, 1) 
# We make the assumption that the number of sons and daughters for the 
# proband and all siblings, is the same. Nieces and nephews of the proband 
# are not sampled separately.
nGrandchild = c(1, 2) 

# Simulate family using `PedUtils` code
fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                    PanelPRODatabase, genes, cancers, 
                    includeGeno = FALSE, includeBiomarkers = TRUE)
                    
# PanelPRO can be run on the simulated family
out = PanelPRO:::PanelPRO11(fam)
```

- The family simulations have several limitations that could be expanded on, including: 
    - Functions simulate a four-generation family with the proband in the third generation. You can make the family smaller by removing relatives, and you can change who the proband is, but there isn't currently a way to specify this in the function. There also isn't a clear-cut way to make a bigger family. 
    - Twins, contralateral breast cancer, and interventions/surgeries are not simulated (the family assumes nobody is a twin, gets contralateral breast cancer, or has an intervention)
    - Age ranges are currently hard-coded, with the first generation being centered at 85 and everybody else simulated based on that. This could be made more flexible/allow user inputs. 
    - You cannot simulate a family in which individuals have different prevalences/penetrances. 
    - Ancestry is automatically set to "nonAJ" and race is automatically set to "All_Races" (regardless of what prevalences/penetrances you feed in)
    - The code assumes that having more than maxMut mutations or having mutations on both alleles is non-admissible (such situations occur more frequently if the allele frequencies are large). Currently, there are no underlying safeguards in the simulation to prevent such scenarios from occurring, so the function's approach is to keep trying to re-simulate the genotype matrix (up to maxTries times) until an admissible genotype has been created for every family member. 
    
    
    
Notes for plotting the pedigree
- dependencies: kinship2, RColorBrewer

```
library(PanelPRO)
plotPed(small.fam)
```

