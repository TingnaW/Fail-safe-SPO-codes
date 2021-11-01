# Fail-safe SPO
Here are the **data and codes** to reproduce the results in the paper ***On improved fail-safe sensor distributions for a structural health monitoring system***.

- A series of tests were performed on a **glider wing** in an environmental chamber to provide a data set suitable for this research, which is in the folder **`ModeShapeData`**. 
- The deployment of 36 candidate sensors is shown in Figure `SensorDeployment.jpg`.

### `1ModeSExtraction`
Codes in this folder can be used to extract features for the next step. The adopted features include the mode shapes and corresponding labels of structural states.
1. **`MSExtraction_ModalIdentification.m`** will provide features for **modal identification**. Two criteria are considered, including the DFIM and the DFIMADPR.
2. **`MSExtraction_DamageDetection.m`** will provide features for **damage detection**. Two criteria are considered, including the SSC and the SSCADPR.

### `2DataProcessing`
1. The codes of **the classical SPO, the fail-safe SPO and the fail-safe with redundancy** SPO corresponding to **four optimisation objectives**, including the DFIM, DFIM-ADPR, SSC and SSCADPR, are in folders **`FIM`, `FIMADPR`, `SSC` and `SSCADPR`** seperately.
   - Two search algorithms were applied, including a deterministic algorithm–exhaustive search (ES) and a stochastic algorithm–genetic algorithm (GA), to obtain the optimal sensor layouts.
   - ExhaustiveS refers to ES and SPOGA refers to GA. 
   - The file with _Improved_ in its name should be run _last_ in the corresponding subfolder to provide the improved fail-safe SPO results.
2. **`FIM-FIMADPR-Assessment`** and **`SSC-SSCADPR-Assessment`** contain performance evaluation codes. For the specific related questions, please refer to Section 5.3 of the paper. 
3. **`OptimalSensorDistribution`** contains codes that facilitate comparisons of optimal sensor distributions. The corresponding part in the paper is Section 5.4.
