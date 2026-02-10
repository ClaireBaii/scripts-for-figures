# scripts-for-figures

Code and data-processing scripts for paper figures (Figure 2 to Figure 9).

## Repository Layout

- `data/`: input datasets used by scripts.
- `figure2/` ... `figure9/`: figure scripts and outputs.
- `environment.yaml`: conda environment definition.
- `setup_env.sh`, `setup_env.ps1`: environment setup and package verification.

## 1) Environment Setup

Run from repository root.

### macOS / Linux

```bash
bash setup_env.sh
```

### Windows (PowerShell)

```powershell
pwsh -File .\setup_env.ps1
```

Both scripts will:

1. Update/create the conda environment defined in `environment.yaml`.
2. Verify required R packages.

## 2) Reproduce Figures

Activate the conda environment first:

```bash
conda activate paper_fig_env
```

Then run scripts from repository root:

```bash
Rscript "figure2/fig2_make_main.R"
Rscript "figure3/figure3_main_alpha_pcoa.R"
Rscript "figure3/figure3_distance_to_control.R"
Rscript "figure3/figure3_rarefaction_robustness.R"
Rscript "figure4/figure4 单样本分组热图.R"
Rscript "figure5/fig5_make_all.R"
Rscript "figure5/fig5_threshold_sensitivity_table.R"
Rscript "figure6/Fig6_phylum_trajectory_with_stats.R"
Rscript "figure7/fig7_segmented_regression.R"
Rscript "figure8/figure8 Procrustes 分析与 PROTEST 检验（美化）.R"
Rscript "figure9/figure9_make_fig9.R"
```

Outputs are written into the corresponding `figure*/` folders.

## Notes

- Temporary diagnostics, local package caches, and run logs are intentionally excluded from version control.
- If your local conda environment name differs, pass it explicitly in PowerShell:

```powershell
pwsh -File .\setup_env.ps1 -EnvName your_env_name
```
