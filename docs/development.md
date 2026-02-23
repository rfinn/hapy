# Architecture Overview

hapy is structured in three layers:

- masktools → low-level masking logic (no GUI, no catalogs)
- hatools → Hα science pipeline (headless, batchable)
- hagui → GUI wrapper around hatools

MaskEngine owns all mask state.
HASession owns science state.
GUI never performs scientific computation directly.
