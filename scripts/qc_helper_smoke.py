import numpy as np
from astropy.table import Table

from qc_helpers import (
    safe_bool_array,
    safe_float_array,
    safe_str_array,
    add_science_columns,
    add_qc_columns,
    add_qc_tier,
    prepare_analysis_table,
)

tab = Table()

tab["VFID"] = ["VFID0001", "VFID0002", None]

tab["PHOT_OK"] = [True, "False", "yes"]
tab["R_PROFILE_OK"] = [True, True, False]
tab["H_PROFILE_OK"] = [True, False, True]
tab["HAPY_MORPH_OK"] = [True, True, False]
tab["H_SM_OK"] = [True, False, True]
tab["GAL_NC_OK"] = [True, False, False]
tab["GAL_CV_OK"] = [False, False, True]

tab["BRIGHT_STAR_FLAG"] = [False, True, False]
tab["ELL0_MASK_WARN"] = [False, False, True]
tab["ELL_MISMATCH"] = [False, False, True]

tab["FILTER_CORRECTION"] = [1.05, 1.35, np.nan]

tab["R50_ARCSEC"] = [12.0, 8.0, np.nan]
tab["H50_ARCSEC"] = [10.0, np.nan, 4.0]
tab["R75_ARCSEC"] = [20.0, 15.0, np.nan]
tab["H75_ARCSEC"] = [16.0, np.nan, 8.0]
tab["R25_ARCSEC"] = [25.0, 18.0, 12.0]
tab["H_MAXDET_ARCSEC"] = [18.0, 4.0, np.nan]

tab["R_PETRO_R50_ARCSEC"] = [11.0, 7.5, np.nan]
tab["H_PETRO_R50_ARCSEC"] = [9.0, np.nan, 3.0]

tab["R_HAPY_GINI"] = [0.60, 0.45, 0.30]
tab["H_HAPY_GINI"] = [0.70, 0.20, np.nan]
tab["R_HAPY_M20"] = [-1.8, -1.4, -1.1]
tab["H_HAPY_M20"] = [-1.5, -2.0, np.nan]
tab["R_HAPY_ASYM"] = [0.10, 0.25, 0.40]
tab["H_HAPY_ASYM"] = [0.18, 0.80, np.nan]

tab["R_PROFILE_NGOOD"] = [30, 25, 5]
tab["H_PROFILE_NGOOD"] = [12, 4, 2]
tab["R_PROFILE_MASKFRAC_MAX"] = [0.1, 0.2, 0.6]
tab["H_PROFILE_MASKFRAC_MAX"] = [0.1, 0.4, 0.7]
tab["H_HAPY_NPIX"] = [200, 20, 5]
tab["H_HAPY_SNP_DET"] = [8.0, 2.0, 1.0]

out = prepare_analysis_table(tab)

print(out.colnames)
print(out["VFID", "DELTA_GINI", "DELTA_M20", "DELTA_ASYM",
          "H50_R50_RATIO", "R_STRUCTURE_GOOD", "HA_EXTENT_GOOD",
          "HA_MORPH_GOOD", "SCIENCE_READY", "QC_TIER"])

print("safe_book_array:",safe_bool_array(tab, "PHOT_OK"))
print("safe_float_array:",safe_float_array(tab, "FILTER_CORRECTION"))
print("safe_str_array:",safe_str_array(tab, "VFID"))

assert "DELTA_GINI" in out.colnames
assert "DELTA_M20" in out.colnames
assert "DELTA_ASYM" in out.colnames
assert "H50_R50_RATIO" in out.colnames

assert np.isclose(out["DELTA_GINI"][0], 0.10, equal_nan=False)
assert np.isclose(out["H50_R50_RATIO"][0], 10.0 / 12.0, equal_nan=False)


for col in [
    "R_STRUCTURE_GOOD", "HA_EXTENT_GOOD", "HA_MORPH_GOOD",
    "SCIENCE_READY", "SCIENCE_PROBLEM", "FILTER_WARNING"
]:
    assert col in out.colnames, f"missing {col}"


print(out["QC_TIER"])



##  TESTING ON REAL TABLE
print("\nNow testing on real table...\n")

tab = Table.read("merged_results.FITS")
out = prepare_analysis_table(tab)

print(len(tab), len(out))
print("QC_TIER" in out.colnames)
print("SCIENCE_READY" in out.colnames)
print("DELTA_GINI" in out.colnames)
print("H50_R50_RATIO" in out.colnames)


import numpy as np

tier = np.array(out["QC_TIER"])
for t in ["A", "B", "C", "D", "F"]:
    print(t, np.sum(tier == t))

print("SCIENCE_READY:", np.sum(out["SCIENCE_READY"]))
print("HA_MORPH_GOOD:", np.sum(out["HA_MORPH_GOOD"]))
print("R_STRUCTURE_GOOD:", np.sum(out["R_STRUCTURE_GOOD"]))
