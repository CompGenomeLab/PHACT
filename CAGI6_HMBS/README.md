

**This folder contains the data and code required to reproduce the results of the CAGI6 HMBS challenge.**

The experimental data is obtained from [MaveDB](https://www.mavedb.org/#/experiments/urn:mavedb:00000108-a) [1]. As mentioned in Zhang et al. (2024) [2], mutations with experimental scores over 1.36 were excluded. The comparisons are made over 5,811 missense mutations.

PHACT scores for 3,023 proteins have been available since **18/10/2022** at the following link: [PHACT data](https://aperta.ulakbim.gov.tr/record/240637). Details on how these scores were computed, including the pseudocode of the PHACT approach, can be found in the [PHACT paper](https://doi.org/10.1093/molbev/msac114) [3].


We participated in the CAGI6 HMBS challenge with a draft version of our approach. In the published version, we considered both position diversity and implemented a scaling approach on the final score. These results were used in the challenge to assess algorithm success. The rank-based results are presented in the table below:

<p align="center">
  <img src="images/Table.png" alt="Alt text" width="500"/>
</p>


As shown in the table, PHACT outperformed both our draft algorithm and PolyPhen-2, as well as the baseline approach.

**References**

1. van Loggerenberg, W., Sowlati-Hashjin, S., Weile, J., Hamilton, R., Chawla, A., Sheykhkarimli, D., ... & Roth, F. P. (2023). Systematically testing human HMBS missense variants to reveal mechanism and pathogenic variation. The American Journal of Human Genetics, 110(10), 1769-1786.
2. 2. Zhang J., Kinch L., Katsonis P. Lichtarge O., Jagota M., Song Y., Sun Y., Shen Y., Kuru N., Dereli O., Adebali O., Alladin M.A., Pal D., Capriotti E., Turina M.P., Savojardo C., Martelli P. L., Babbi G., Bovo S., Hu Z., Pucci F., Rooman M., Cia G., Tsishyn M., Strokach A., van Loggerenberg W., Roth F.P., Radivojac P., Brenner S.E., Cong Q., Grishin N.V. (2024). Assessing predictions on fitness effects of missense variants in HMBS in CAGI6.
3. Kuru N., Dereli O., Akkoyun E., Bircan A., Tastan O., & Adebali O. (2022). PHACT: Phylogeny-aware computing of tolerance for missense mutations. Molecular Biology and Evolution, 39(6), msac114.





