# Online Supplement to Mansfield & Avestruz (2020)

This repository contains supplementary data for the paper Mansfield & Avestruz (2020). This includes catalogues of convergence limits at various tolerances (delta = 0 to 0.10), polynomial fits to the z=0 mean mass relations, <X(M_vir)> for various halo properties, and code which implements the bias model implemented in Section 5.

`convergence_limits/` contains convergence limits for several simulations. This is the M_vir at which this simulaltion deviates from a mean relation to within a given fracitonal accuracy, delta. To find the delta=0.05 limits for subhaloes for the variable V_max, use the file `convergence_limits/vmax/sub/vmax_sub_delta05.txt`. The contents of these files are described in a header comment.

`fits/` contains polynomial fits to various mass relations. 

`scripts/` contains scripts which implement the mass profile debiasing model described in Section 5 of the base paper.
