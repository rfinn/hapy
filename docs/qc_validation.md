# Architecture

```
qc_helpers.py
    ↓ shared logic

qc_results.py
    → applies qc_helpers
    → outputs flags, outliers, QC tables

validate_duplicates.py
    → uses duplicate helpers + stats

validate_measurements.py
    → uses derived columns + stats

validate_dashboards.py
    → reads outputs, builds visuals

science_firstlook.py
    → uses cleaned/filtered data only
```

# Ownership

```
qc_results.py            -> QC triage, flags, suspicious objects
validate_duplicates.py   -> repeatability from duplicate observations
validate_measurements.py -> reliability/stability of metrics more generally
validate_dashboards.py   -> dashboards/review products
science_firstlook.py     -> science plots using trusted subsets
```
