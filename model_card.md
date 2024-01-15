# Model Card: HSMA4-12-DES-rheum

## Model Details

The implementation of HSMA4-12-DES-rheum-mkdocs model within this repository was created as part of a Health Service Modelling Associates 4 (HSMA4) project undertaken by Martina Fonseca and Xiaochen Ge. This model card describes the updated version of the model, released {10-01-2024}, [cb503ec](https://github.com/nhsx/HSMA4-12-DES-rheum/commit/cb503ec9dded1a173334ee1e5cf023f015b28827).

## Model Use

### Intended Use

The model and code is experimental.

This model is intended to be used to generate evidence for understanding the dynamics in the outpatient setting and waiting list management - and the effect of interventions like Advice & Guidance or PIFU. The intended use requires the model inputs and parameters to be calibrated for a specific trust and the results to be validated alongside expert opinion.

Further guidance can be found in the README or the (internal) MkDocs documentation.

### Out-of-Scope Use Cases

This model is a simple representation of the outpatient and waiting list setting.   Currently no dynamic escalations or human behaviour characteristics are included in the model.  Without these formal and informal mitigation dynamics, the model is overly sensitive and shows an explosive nature when resources become depleted.   Therefore, the results from this model should be used as an indication to inform the narrative and not used for performance monitoring or to be used in isolation for planning purposes.
