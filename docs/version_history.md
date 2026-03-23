## Meaning of Version Increments

#### Patch-level examples (+0.01)
- fix a bad filename
- fix one flag not being written correctly
- improve logging
- rerun after a narrow implementation bug

#### Minor-version examples (+0.10)
- fix units on statmorph radii
- add derived columns
- change how profiles are measured
- add corrected radii
- change masking behavior that affects science outputs

### Major-version examples (+1.00)
- redesign results table
- replace major parts of morphology/profile pipeline
- break backward compatibility

## v0.2 (2026-03-23)
- fixed pixel vs arcsec inconsistency in radial measurements; column names for radial sizes now explicitly have _ARCSEC
- added vr and redshift to metadata in `get_cutouts` and propagated through to `run_analysis`
- validated telescope-dependent offsets as a unit-conversion issue
- updated QC + validation framework
- flagged need for inclination correction (not yet applied)
