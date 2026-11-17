# Paper II: NPB reproducibility deposit (17 November 2026)

Companion to **Terminal scalar states, Hadamard obstruction and local stress in
the Schwarzschild interior**, Francisco R. Villatoro.

`Paper_II_NPB_supplementary_code.zip` is the exact frozen NPB code package.
It contains Python sources, historical reference outputs, the independent audit,
and the 17 November 2026 validation report. Its internal `SHA256SUMS` covers
61 files. The outer checksum here identifies the exact ZIP bytes.

This release packages existing calculations and their validation; it introduces
no new scientific algorithms. Historical AoP filenames are retained for provenance.
The assessment inside the ZIP distinguishes numerical residuals from rigorous
bounds. Consult its README for dependencies and reproduction commands.

```bash
sha256sum -c SHA256SUMS
unzip Paper_II_NPB_supplementary_code.zip
cd Paper_II_NPB_supplementary_code
sha256sum -c SHA256SUMS
python3 -m pip install -r requirements.txt
python3 run_npb_checks.py
```

The scientific manuscript and editorial correspondence are not part of this
software deposit. Cite release v1.3.1 and its version-specific Zenodo DOI once
assigned; older version DOIs do not identify this expanded NPB package.
