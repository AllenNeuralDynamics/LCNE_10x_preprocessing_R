"""
Add AIND metadata for the ``LCNE-transcriptomics-preprocessing-from-R`` data asset.

This asset is scientist-derived data: Allen-internal locus-coeruleus norepinephrine
(LC-NE) transcriptomics that was preprocessed in **R** (Seurat) and imported into the
``LCNE_transcriptomics_preprocess`` analysis capsule. It contains:

* ``snRNAseq_LCNE.h5ad`` -- single-nucleus RNA-seq (10xV4, mouse C57BL6J, ROI
  "Mouse 10x PONS - LC", study ``Neuromodulatory_Noradrenergic``, multiple donors,
  M/F batches).
* ``retroseqdata_raw_from_R/`` -- retro-seq raw count matrix (``retro_raw.mtx``,
  ``cells.tsv``, ``genes.tsv``, ``metadata.csv``).

Unlike the external Nardone 2024 MERFISH asset (see ``add_metadata.py``), this is
AIND-generated derived data, so we build a *standalone* AIND ``DataDescription``
(institution = AIND, data_level = derived) plus a ``Processing`` object documenting the
R-preprocessing step that produced it. Because it aggregates many donors across two
modalities with no single AIND base asset, no base metadata is inherited.

Every provenance identifier is **hard-coded** below so this script is self-contained and
can be run from anywhere -- it does NOT read ``.codeocean/datasets.json`` or local git,
which would not be available outside the original capsule. The values are:

* the **data asset** mount name + Code Ocean asset id (of this asset), and
* the **code capsule / commit hash** of ``LCNE_10x_preprocessing_R`` -- the separate R
  capsule that actually produced the data -- pinned to the commit that generated it.

If any of these change, edit the constants in the CONFIG block.

Run::

    python add_metadata_LCNE.py [output_dir]

This writes only ``data_description.json`` and ``processing.json`` (no data copy); by
default to ``/results/LCNE-transcriptomics-preprocessing-from-R``. Attach the asset's data
alongside these JSONs when creating the Code Ocean data asset, then file an issue to
transfer it to ``aind-open-data`` before publication.
"""

import sys
from datetime import datetime
from pathlib import Path

import aind_data_schema.core.data_description as ds
import aind_data_schema.core.processing as ps
from aind_data_schema.core.metadata import Metadata

# --------------------------------------------------------------------------------------
# CONFIG -- all identifiers hard-coded so the script is portable (no capsule/git lookups).
# --------------------------------------------------------------------------------------

# The data asset this metadata describes.
MOUNT_NAME = "LCNE-transcriptomics-preprocessing-from-R_2026-07-08_11-11-11"
ASSET_ID = "9b928b21-b5ae-4cab-90e7-ec68e78a327c"
CREATION_TIME = datetime(2026, 7, 8, 11, 11, 11)  # matches the mount timestamp suffix

# The code that produced the asset: the *separate* R preprocessing capsule
# ``LCNE_10x_preprocessing_R`` (NOT this analysis capsule). COMMIT_HASH is the commit that
# generated the asset -- the R repo's state on 2026-07-08, before the asset was created.
CODE_NAME = "LCNE_10x_preprocessing_R"
REPO_URL = "https://github.com/AllenNeuralDynamics/LCNE_10x_preprocessing_R"
COMMIT_HASH = "bdf73dfa38470072d672004ba75ba704365df4e4"  # 2026-07-08, "Edited 03_LC_subclustering.R"

CO_ASSET_URL_BASE = "https://codeocean.allenneuraldynamics.org/data-assets"

# --------------------------------------------------------------------------------------


def build_processing() -> ps.Processing:
    """Document the R-preprocessing step (LCNE_10x_preprocessing_R) that produced this asset."""
    run_time = CREATION_TIME.strftime("%Y-%m-%dT%H:%M:%S")
    code_details = ps.Code(
        name=CODE_NAME,
        url=REPO_URL,
        version="1.0",
        commit_hash=COMMIT_HASH,
        input_data=[
            ps.DataAsset(name=MOUNT_NAME, url=f"{CO_ASSET_URL_BASE}/{ASSET_ID}"),
        ],
    )
    return ps.Processing(
        data_processes=[
            ps.DataProcess(
                process_type=ps.ProcessName.OTHER,
                name="transcriptomics_preprocessing_from_R",
                stage=ps.ProcessStage.PROCESSING,
                experimenters=["Shuonan Chen"],
                start_date_time=run_time,
                end_date_time=run_time,
                code=code_details,
                notes=(
                    "LC-NE snRNA-seq and retro-seq were preprocessed in R (Seurat) in the "
                    "LCNE_10x_preprocessing_R capsule: raw count matrices assembled and "
                    "exported as MTX/TSV and AnnData, with quality control, normalization, "
                    "and filtering to the LC-NE population. The resulting datasets were then "
                    "imported into the LCNE_transcriptomics_preprocess analysis capsule for "
                    "downstream analysis. The recorded url/commit_hash point to "
                    "LCNE_10x_preprocessing_R at the commit that generated this asset."
                ),
            ),
        ],
    )


def build_data_description() -> ds.DataDescription:
    """Standalone AIND derived-data description for the LC-NE transcriptomics asset."""
    return ds.DataDescription(
        name=MOUNT_NAME,
        creation_time=CREATION_TIME,
        institution=ds.Organization.AIND,
        data_level=ds.DataLevel.DERIVED,
        investigators=[ds.Person(name="Shuonan Chen")],
        project_name="LC-NE transcriptomics",
        # scRNAseq is the closest available modality; there is no dedicated
        # single-nucleus RNA-seq / retro-seq modality in the schema registry.
        modalities=[ds.Modality.SCRNASEQ],
        license=ds.License.CC_BY_40,
        # The schema requires >=1 funding_source. Grant number omitted (matches the
        # Nardone script); AIND is recorded as the funder.
        funding_source=[ds.Funding(funder=ds.Organization.AIND)],
        data_summary=(
            "Allen Institute for Neural Dynamics locus-coeruleus norepinephrine (LC-NE) "
            "transcriptomics (study Neuromodulatory_Noradrenergic, mouse C57BL6J, ROI "
            "'Mouse 10x PONS - LC'). Contains single-nucleus RNA-seq (10xV4; "
            "snRNAseq_LCNE.h5ad) and retro-seq raw count matrices "
            "(retroseqdata_raw_from_R/: retro_raw.mtx, cells.tsv, genes.tsv, metadata.csv) "
            "across multiple donors and M/F batches. Preprocessed in R (Seurat): raw "
            "count assembly, quality control, normalization, filtering, and MTX/TSV -> "
            "AnnData export, then imported into the LCNE_transcriptomics_preprocess "
            "capsule for downstream analysis."
        ),
    )


DEFAULT_OUTPUT = "/results/LCNE-transcriptomics-preprocessing-from-R"


def main(output_path: str = DEFAULT_OUTPUT) -> None:
    Path(output_path).mkdir(parents=True, exist_ok=True)

    metadata = Metadata(
        name=MOUNT_NAME,
        location=f"s3://aind-open-data/{MOUNT_NAME}",
        data_description=build_data_description(),
        processing=build_processing(),
    )
    metadata.data_description.write_standard_file(output_path)
    metadata.processing.write_standard_file(output_path)

    print(f"Asset mount  : {MOUNT_NAME}")
    print(f"Asset id     : {ASSET_ID}")
    print(f"Code capsule : {CODE_NAME} ({REPO_URL})")
    print(f"Commit hash  : {COMMIT_HASH}")
    print(f"Wrote data_description.json and processing.json to {output_path}")


if __name__ == "__main__":
    positional = [a for a in sys.argv[1:] if not a.startswith("-")]
    out = positional[0] if positional else DEFAULT_OUTPUT
    main(out)
