# Uploading the code package to GitHub

Repository: <https://github.com/Alexios130/mannalen>

Archive to extract: `WorstScenarioPath_Code.zip`

## Upload through the GitHub website

1. Download and extract `WorstScenarioPath_Code.zip` on your computer.
2. Open the extracted folder and upload its **contents** to the repository's
   `main` branch. Do not upload only the ZIP: reviewers must be able to browse
   `README.md`, `certificate/`, `experiments/`, and `verification/` directly.
3. Confirm in a signed-out browser that at least these paths are public:

   - `README.md`;
   - `REPRODUCIBILITY.md`;
   - `certificate/verify_n3_certificate.py`;
   - `certificate/n3_certificate.json`;
   - `experiments/REPRODUCIBILITY_MANIFEST.md`;
   - `verification/arb_verification_output.txt`.

4. Confirm that all nine disclosed Python programs listed below are present:

   - `certificate/verify_n3_certificate.py`;
   - `certificate/generate_n3_certificate.py`;
   - `certificate/independent/rational_interval_check.py`;
   - `certificate/independent/compute_decimal.py`;
   - `experiments/search_switches.py`;
   - `experiments/linear_family_search.py`;
   - `experiments/check_linear_family_decimal.py`;
   - `experiments/grid_resolution_study.py`;
   - `experiments/make_n3_path_figure.py`.

5. Confirm that the GitHub upload contains **none** of the following:

   - the manuscript PDF;
   - LaTeX source or bibliography files;
   - literature-audit or internal-review files;
   - old repository, supplementary, or legacy ZIP archives.

6. Run the commands in `README.md` once from a clean checkout.
7. Create a tag or GitHub release, for example `v1.0.0`, for the exact version
   supplied to the reviewers, and keep that URL with the submission records.

## Files submitted separately to the journal

The compiled manuscript PDF and editable LaTeX source ZIP are uploaded through
the journal portal. They do not belong in `WorstScenarioPath_Code.zip` and do
not need to be placed in the public code repository.

## Licence and permanent archive

No software licence has been selected. Without one, the public files remain
under ordinary copyright. A GitHub tag gives the reviewed code version a fixed
name. Before final publication, that tag can also be archived in a service such
as Zenodo and assigned a DOI, so later changes to `main` do not alter the cited
artifact. A DOI is not required to begin peer review.
