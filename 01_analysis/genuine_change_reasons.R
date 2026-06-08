
# Genuine changes - tables and upset chart for comparison

library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(eulerr)
library(ggthemes)
library(ggupset)
library(tidyr)

work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"




# summarise primary and secondary threats from RLI assessments

# no way to pull these out so has to be manual summary
# 12 species

# 1. Acianthera hondurensis (Ames) Pridgeon & M.W.Chase
# Reduction of area of habitat likely >10%
# cattle ranching, urban development

# 2. Aeschynomene kerstingii Harms
# increases in agricultural activities combined with human settlement
# and population growth.

# 3. Aloe pubescens
# Habitat encroachment by expanding agriculture is the main threat
# type location

# 4. Astrocaryum carnosum F.Kahn & B.MillÃ¡n
# The 'Palmas del Espino' palm oil plantation is likely to have affected
# the population of this species along with ongoing agricultural expansion
# such a coca production, and logging, along the river valley

# 5. Bulbophyllum multiflorum
# the mining project in Ambatovy, which is the main threat to the species
# has caused the status to deteriorate enough to move it into a higher threat category (VU)
# Driver of change or threat mitigated *
# Primary: Energy production and mining
# Most Important Secondary: Logging and wood harvesting

# 6. Dalbergia glomerata
# Since first assessment of the species, there has been further population decline
# occurring at a high rate due to an increase in harvesting of rosewood species.
# This is driven by greater international demand for this timber type, leading to
# legal and illegal logging activities which target Dalbergia spp.. The range of the species is now also considered smaller as only Mexico is now included in it's native range.

# 7. Dalbergia peltieri
# The exploitation of Dalbergia has increased in the last ten years

# 8. Neltuma kuntzei (Harms) C.E.Hughes & G.P.Lewis
#Logging and conversion of land for agriculture

# 9. Neoharmsia madagascariensis
#Loss of existing habitat due to increased pressure slash and burn, shifting agriculture

# 10. Newtonia erlangeri (Harms) Brenan
# Charcoal burning, logging, habitat loss

# 11. Xyris boliviana
#mining upstream

#####################################################
# 4. Dypsis scottiana (Becc.) Beentje & J.Dransf.
# removed, no longer changing. trumped by criterion B
#####################################################




