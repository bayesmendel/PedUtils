# PedUtils
simulations+model evaluations

Dependencies: `PanelPRO`, `abind`

Notes from Jane 6/19/2020
- I added my family simulation functions this repo (prefixed with `sim.`). They are compatible with (and dependent on) the current PanelPRO package. 
- Since the PanelPRO package is still being updated, the simulation functions may also need updates. Please let me know if you run into any issues. 
- Example usage (need to first run PanelPRO's `buildDataBase` and `calc_CP` functions on a dummy family to obtain the densty and survival functions for the appropriate cancers/genes): 

```
library(PanelPRO)

## 11 cancers and 11 genes to use
# Short cancer names (including CBC)
cancers = c("BRA", "BC", "COL", "ENDO", "GAS", "KID", 
            "MELA", "OC", "PANC", "PROS", "SMA", "CBC")
# Look up long cancer names (don't include CBC)
cancers_long = PanelPRO:::CANCER_NAME_MAP$long[sapply(cancers[-length(cancers)], function(x){
  which(x==PanelPRO:::CANCER_NAME_MAP$short)
})]
# Genes
genes = c("ATM", "BRCA1", "BRCA2", "CDKN2A", "CHEK2", "EPCAM", 
          "MLH1", "MSH2", "MSH6", "PALB2", "PMS2")


## Create a dummy family to generate penetrance densities and survivals
# Empty data frame of cancer affectation statuses
dummy.cancers = setNames(as.data.frame(matrix(0, nrow=2, ncol=length(cancers))), 
                         paste0("isAff", cancers))
# Empty data frame of cancer ages
dummy.ages = setNames(as.data.frame(matrix(NA, nrow=2, ncol=length(cancers))), 
                      paste0("Age", cancers))
# Dummy family
dummy.fam = data.frame(ID=c(1,2), 
                       MotherID=c(0,1), 
                       FatherID=c(0,1), 
                       Sex=c(0,1), 
                       isProband=c(0,1), 
                       Twins=c(0,0), 
                       Ancestry=rep("nonAJ", 2), 
                       CurAge=c(60,30), 
                       isDead=c(0,0), 
                       race=rep("All_Races", 2), 
                       dummy.cancers, 
                       dummy.ages)
# Assign no risk modifiers
dummy.fam$interAge = dummy.fam$riskmod = list(character(0))


## Get cancer penetrance densities and survivals from the dummy family
# Build a dummy database
dummy.db = buildDatabase(genes=genes, 
                         cancers=cancers_long)
dummy.db$Contralateral = PanelPRODatabase$Contralateral

# Run `checkFam` on the dummy family
dummy.fam.checked = checkFam(dummy.fam, dummy.db)$ped_list[[1]]

# Cancer penetrance densities and survivals
CP = calcCancerPenetrance(dummy.fam.checked, dummy.db, 
                          max_mut=2, net=TRUE, consider.modification=FALSE)

# Extract allele frequencies from database
alleleFreq = PanelPRODatabase$AlleleFrequency[,"nonAJ"][genes]


## Simulate family

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
fam = sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                 alleleFreq, CP, genes, cancers_long, includeGeno=FALSE)
# PanelPRO can be run on the simulated family
out = PanelPRO:::PanelPRO11(fam)
```

- The family simulations have several limitations that could be expanded on, including: 
    - Functions simulate a four-generation family with the proband in the third generation. You can make the family smaller by removing relatives, and you can change who the proband is, but there isn't currently a way to specify this in the function. There also isn't a clear-cut way to make a bigger family. 
    - Twins, contralateral breast cancer, and risk modifiers are not simulated (the family assumes nobody is a twin, gets contralateral breast cancer, or has a risk modifier/intervention)
    - Age ranges are currently hard-coded, with the first generation being centered at 85 and everybody else simulated based on that. This could be made more flexible/allow user inputs. 
    - You cannot simulate a family in which individuals have different prevalences/penetrances. 
    - Ancestry is automatically set to "nonAJ" and race is automatically set to "All_Races" (regardless of what prevalences/penetrances you feed in)
    
    
    
Notes for plotting the pedigree
- dependencies: kinship2, RColorBrewer

```
library(PanelPRO)
plotPed(small.fam)
```

