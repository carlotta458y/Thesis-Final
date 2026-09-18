# Cost-effectiveness model — specific AI-based medical device

R implementation of the health-economic model used in the MSc thesis: Economic evaluation of AI- based medical devices (Management of Technology, TU Delft, 2026).

The model evaluates **a specific AI-based medical device**, which is an AI-based preoperative planning service for
Endovascular Aneurysm Repair (EVAR) developed by a Medtech company. The tool assigns each
patient an endoleak-risk indication (ERI) before surgery; a high indication can move
the patient to a different surgery type. The model asks whether that re-allocation is
cost-effective over a lifetime horizon, and what price the tool could command.

---

## What the model does

A **Markov cohort model** over seven health states:

`sac_shrinking`, `sac_stable`, `sac_expansion`, `Rupture`, `PostLeak`, `OSR`, `Dead`

- **Perspective** — hospital payer
- **Horizon** — lifetime (25 years)
- **Cycles** — 30 days, then annual to the end of the horizon
- **Half-cycle correction** — QALYs and state costs use the average of start- and
  end-of-cycle occupancy, discounted at mid-cycle
- **Surgery types** — Standard EVAR, FEVAR, Endoanchor, open surgical repair (OSR)
- **Outcome** — Net Monetary Benefit at three willingness-to-pay thresholds
  (€20k, €50k, €80k per QALY)

The AI-based medical device acts on the treatment decision, not on a health state, so its effect is
applied **before** the Markov model rather than inside it. This is organised in three
layers:

| Layer | Question | Where in the code |
|---|---|---|
| 1 | Who *is* a would-leaker, before any tool? | `wouldleaker_share` in the scenario loop |
| 2 | What changes for one patient type moved between the two surgeries of a decision pair? | the four `run_EVAR` cells |
| 3 | How well does the tool classify, and what does that cost? | sensitivity, specificity, compliance, price |

`Markov_Model.R` folds layers 2 and 3 together and reports one cohort result per
scenario. `Accuracy_Model.R` keeps them apart, so the accuracy thresholds and the
headroom price can be solved in closed form. The two routes are algebraically
equivalent; the validation checks in `Accuracy_Model.R` assert this directly.

---

## Repository contents

| File | What it is |
|---|---|
| `Markov_Model.R` | The engine: workbook loading, parameter construction, transition matrices, cohort traces, cost and QALY accumulation, and the base-case scenario loop |
| `Accuracy_Model.R` | Net benefit per patient type, accuracy verdicts and thresholds, value of sensitivity/specificity, value of compliance, headroom price. Includes the validation checks |
| `Confirming_Model.R` | Re-runs the analysis assuming an initial allocation of patients across surgery types, rather than a single source surgery |
| `DSA_Markov.R` | Deterministic sensitivity analysis; tornado diagrams |
| `PSA_MArkov.R` | Probabilistic sensitivity analysis; cost-effectiveness planes |
| `Plot_Markov.R` | Figures: ROC-space cost-effective regions, value of accuracy, value of compliance, headroom, cost-effectiveness plane |
| `parameters_v3.xlsx` | All model inputs (see below) |

---

## Running it

**Requirements** — R 4.5.2, packages `readxl` (1.4.5) and `writexl` (1.5.4).
No other dependencies.

```r
install.packages(c("readxl", "writexl"))
```

Open the repository folder as your working directory, then run whichever analysis you
need. Each script sources the ones it depends on, so running one further down the
chain re-runs everything above it:

```
Markov_Model.R                     base case
  └── Accuracy_Model.R             accuracy, compliance, headroom
        └── Confirming_Model.R     initial-allocation variant
              └── Plot_Markov.R    all figures

Markov_Model.R
  ├── DSA_Markov.R                 deterministic sensitivity analysis
  └── PSA_MArkov.R                 probabilistic sensitivity analysis
```

```r
source("Plot_Markov.R")   # base case + accuracy analysis + every figure
```

Sourcing a file as a dependency defines its functions without re-running its analysis,
so nothing is computed twice.

---

## Inputs

All inputs live in `parameters.xlsx`, one sheet per group:

| Sheet | Contents |
|---|---|
| `general` | Horizon, start age, discount rates, WTP thresholds, tool accuracy, compliance, price |
| `follow_up` | Cumulative event probabilities per surgery type and follow-up moment |
| `life_table` | Annual background mortality by age |
| `scenarios` | The decision pairs, their direction, and the elicited would-leaker share |
| `costs` | Unit costs, each with source and price year |
| `utilities` | Health-state utilities |
| `time` | Procedure, ICU and hospital durations; follow-up intensities |
| `event_probs` | One-off event probabilities |
| `sensitivity_analysis` | Ranges and standard errors driving the DSA and PSA |

Every cost row carries its source and year. Costs are drawn from the Dutch costing
manual, NZa tariffs and published literature.

---

## Outputs

| File | Produced by |
|---|---|
| `results.xlsx` | `Markov_Model.R` — incremental cost, QALYs, endoleaks avoided, NMB, headroom per scenario × threshold |
| `results_temp.xlsx` | `Accuracy_Model.R` — the above plus verdicts, accuracy and compliance thresholds, headroom under three settings |
| `results_confirming.xlsx` | `Confirming_Model.R` |
| `DSA_tornado_*k.pdf` | `DSA_Markov.R` |
| `PSA_ce_plane.pdf` | `PSA_MArkov.R` |
| `plots/*.png` | `Plot_Markov.R` |




## Contact

c.lichtenauer@student.tudelft.nl
