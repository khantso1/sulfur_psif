# sulfur_psif

Contents: epsif_calc.m

Code for a steady-state equilibrium model of bulk and position-specific sulfur isotope fractionation in the S8 (elemental sulfur)-sulfide-polysulfide-pyrite pool. This script calculates a net S8-pyrite isotopic offset using input values q and n. It assumes that pyrite forms by the polysulfide reaction detailed by Luther (1991) and Butler et al. (2004). PSIF calculations for polysulfide chain lengths S4 through S7 use d34S values from Amrani et al. (2006), who measured the d34S of polsyulfide chains sorted by length and obtained distinct d34S values for S4, S5, S6, and S7. PSIF calcs for S9 and its disproportionation to S8 use the 3.4 to 5.4 per mil fractionation observed by Amrani and Aizenshtat (2004).

epsif_calc.m was used to generate Figure 6 in Hantsoo et al. (2023). A copy of epsif_calc.m is also available under the name position_specific_sulfur_isotope_model.m at the Johns Hopkins Research Data Repository (https://doi.org/10.7281/T1/EHSEV2).


References:

Amrani, A., & Aizenshtat, Z. (2004). Mechanisms of sulfur introduction chemically controlled: d34S imprint. Organic Geochemistry, 35(11-12), 1319-1336. https://doi.org/10.1016/j.orggeochem.2004.06.019

Amrani, A., Kamyshny, A., Lev, O., & Aizenshtat, Z. (2006). Sulfur stable isotope distribution of polysulfide anions in an (NH4)2Sn aqueous solution. Inorganic Chemistry, 45(4), 1427-1429. https://doi.org/10.1021/ic051748r
 
Butler, I. B., Böttcher, M. E., Rickard, D., & Oldroyd, A. (2004). Sulfur isotope partitioning during experimental formation of pyrite via the polysulfide and hydrogen sulfide pathways: implications for the interpretation of sedimentary and hydrothermal pyrite isotope records. Earth and Planetary Science Letters, 228(3-4), 495-509. https://doi.org/10.1016/j.epsl.2004.10.005

Hantsoo, K., Gomes, M., Malkin, S., Brenner, D., & Kenney, W. F. (2023). Sedimentary Pyrite Formation in a Seasonally Oxygen‐Stressed Estuary: Potential Imprints of Microbial Ecology and Position‐Specific Isotope Fractionation. Journal of Geophysical Research: Biogeosciences, 128(4), e2022JG007324. https://doi.org/10.1029/2022JG007324.

Luther III, G. W. (1991). Pyrite synthesis via polysulfide compounds. Geochimica et Cosmochimica Acta, 55(10), 2839-2849. https://doi.org/10.1016/0016-7037(91)90449-f
