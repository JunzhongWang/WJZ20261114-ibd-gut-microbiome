# Data and code for: Between-subject divergence without faster within-subject turnover but with frozen within-host strain drift in the gut microbiome of inflammatory bowel disease — a matched-interval analysis in two longitudinal cohorts

Analysis code, data tables and result files for the manuscript (Microbiome submission).

**Version:** v1.0.1 · **Released:** 2026-11-14
**Code:** Python 3 (standard library; one analysis additionally uses NumPy) + three small R helpers + a GT-Pro cloud runbook · **License:** MIT (code), CC BY 4.0 (derived tables and data matrices)

## 1. What this archive contains

| Path | What it is |
|---|---|
| `pipeline/` | Analysis pipeline for every figure and table in the paper. `ibd_build_cohort.py` builds the matched-interval cohort from the public inputs; the remaining `ibd_*.py` scripts compute R (within/between-subject Bray–Curtis ratio), functional-layer results, strain-marker results, cross-layer contrasts and the SNV summaries; the `*_fig*.py` scripts draw the manuscript figures. `fetch_cmd3_direct.py` / `fetch_pathways.py` / `rda_extract*.R` reproduce the curatedMetagenomicData retrieval chain that produced `data/`. |
| `gtpro/` | Runbook for the strain-level (SNV) extension: GT-Pro genotyping of 1,886 metagenomes against the 53-million-variant catalogue (881 species, 20190723 release) on AWS Spot (see `gtpro/README.md`, `DEPLOY_aws.md`). Includes manifest design, filtering, aggregation, the self-test harness, cloud operation scripts (`aws/`) and the catalogue metadata tables. The 32 GB of raw GT-Pro output is **not** redistributed; `results/ibd2/snv*.csv` hold the derived summaries. |
| `results/ibd2/` | Every numerical result cited in the manuscript: cohort tables, R per group and per interval band (`r_*.csv`, `marker_*.csv`, `functional_*.csv`), robustness checks (`robustness.csv`), SNV layer results (`snv*.csv`), flare/medication contrasts, cross-layer table, all figure sources (`figs/`, `png/`) and per-subject marker presence matrices (`marker/`). `README_snv.md` documents the SNV result files. |
| `data/HMP_2019_ibdmdb/` | **Input data** for the Multi'omics cohort (PRJNA398089, cMD3 reprocessing): `species.csv`/`species.tsv` (samples × species, comma/tab), `meta.csv`, `pathway.tsv` (HUMAnN 3 pathway abundance) and the unstratified pathway extract. |
| `data/HallAB_2017/` | Same files for the validation cohort (PRJNA385949). |
| `data/cache/` | Retrieval caches needed by two scripts: `sample_metadata.csv` (curatedMetagenomicData sample metadata) and `mp_*.rda` (MetaPhlAn marker-presence resources used by `export_marker_presence.R`). |
| `MANIFEST.csv` | Every file with byte size and SHA-256, for integrity checking. |
| `CITATION.cff` | Machine-readable citation metadata. |

## 2. Cohorts

| Cohort | BioProject | Design | Samples | Subjects |
|---|---|---|---|---|
| Multi'omics (IBDMDB) | PRJNA398089 | stool every 1–2 weeks, up to 1 year | 1,627 | 130 (CD / UC / non-IBD) |
| Hall 2017 validation | PRJNA385949 | stool monthly, up to 12 months | 259 | 33 (IBD / healthy) |

Abundance and pathway profiles are curatedMetagenomicData v3 reprocessings (bioBakery 3 / MetaPhlAn 3 / HUMAnN 3) of the public cohorts. Raw reads are not redistributed; they remain under the BioProject accessions above. The Multi'omics clinical metadata release (`hmp2_metadata_2018-08-20.csv`) is distributed through https://ibdmdb.org. No material is under embargo.

## 3. What is *not* included

- **Raw GT-Pro output (32 GB, `gtpro/aws_out/`, `gtpro/gtpro_out/`)** — per-SNV genotype calls for 1,886 metagenomes. The derived summaries in `results/ibd2/snv*.csv` are included; `gtpro/README.md` and `DEPLOY_aws.md` document how to regenerate them from the SRA runs of the two BioProjects using the GT-Pro catalogue.
- **curatedMetagenomicData pathway `.rda` caches** — `pipeline/fetch_pathways.py` re-downloads them from ExperimentHub; `data/<study>/pathway.tsv` is included directly, so this step is only needed to rebuild from source.

## 4. Run order

```bash
# ---- stage 1: cohort + core R (minutes) --------------------------------------
python3 pipeline/ibd_build_cohort.py         # -> results/ibd2/cohort_*.csv, abund_*.csv
python3 pipeline/ibd_strata.py               # -> r_by_group.csv, r_pooled.csv
python3 pipeline/ibd_secondary.py            # -> diversity.csv, robustness.csv, subject_*.csv
python3 pipeline/ibd_crosslayer.py           # -> crosslayer.csv (+ summary txt)

# ---- stage 2: functional layer (HUMAnN 3 pathways) ---------------------------
python3 pipeline/ibd_functional.py --boot 500   # -> functional_*.csv
python3 pipeline/ibd_pathways.py                # -> pathway_targeted.csv, pathway_variability.csv
python3 pipeline/ibd_flare_meta.py              # -> flare_meta.csv

# ---- stage 3: strain markers (MetaPhlAn marker presence) ----------------------
Rscript pipeline/export_marker_presence.R HMP_2019_ibdmdb <taxid list>   # -> marker/*.csv.gz
Rscript pipeline/export_marker_presence.R HallAB_2017 <taxid list>
python3 pipeline/ibd_strain_markers.py          # -> marker_R*.csv, marker_by_band*.csv
python3 pipeline/ibd_pairs_analyses.py          # -> subject_tests.csv, turnover_slopes.csv  [needs numpy]

# ---- stage 4: SNV layer (requires external GT-Pro run, see gtpro/README) ------
python3 pipeline/ibd_snv_collect.py             # -> snv*.csv        [reads GT-Pro aggregation]
python3 pipeline/ibd_snv_summary_fig.py         # -> fig_snv_summary.svg   (= figs/fig7_strain.svg)
python3 pipeline/ibd_snv_persistence_fig.py     # -> fig_snv_persistence.svg (= figs/fig8_persistence.svg)
python3 pipeline/ibd_snv_robustness_fig.py      # -> fig_snv_robustness.svg (= figs/fig9_strain_robust.svg)

# ---- stage 5: figures (< 1 min) ------------------------------------------------
python3 pipeline/ibd_figs.py     # -> figs/fig1_design.svg ... figs/fig5_robust.svg
python3 pipeline/ibd_fig6.py     # -> figs/fig6_function.svg
```

Figure → script mapping: Fig 1–5 `ibd_figs.py`; Fig 6 `ibd_fig6.py`; Fig 7 `ibd_snv_summary_fig.py`; Fig 8 `ibd_snv_persistence_fig.py`; Fig 9 `ibd_snv_robustness_fig.py`.

## 5. Licence

Code: MIT. Derived tables and data matrices (`data/`, `results/`): CC BY 4.0. If you reuse the tables, please cite the article and this archive.

## 6. Contact

Junzhong Wang, Department of Infectious Diseases, Union Hospital, Tongji Medical College, Huazhong University of Science and Technology — wangjunzhong@hust.edu.cn
