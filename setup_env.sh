#!/usr/bin/env bash
set -euo pipefail

ENV_FILE="${1:-environment.yaml}"

if [[ ! -f "$ENV_FILE" ]]; then
  echo "Environment file not found: $ENV_FILE" >&2
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found in PATH. Install Miniforge/Miniconda and run 'conda init' first." >&2
  exit 1
fi

ENV_NAME="$(awk -F': *' '/^name:/{print $2; exit}' "$ENV_FILE")"
if [[ -z "${ENV_NAME:-}" ]]; then
  echo "Failed to parse environment name from $ENV_FILE" >&2
  exit 1
fi

echo "Updating conda environment '$ENV_NAME' from '$ENV_FILE'..."
conda env update -f "$ENV_FILE" --prune

echo "----------------------------------------------------------------"
echo "Verifying required R packages in '$ENV_NAME'..."
echo "----------------------------------------------------------------"

verify_script="$(mktemp)"
cat >"$verify_script" <<'RSCRIPT'
pkgs <- c(
  "ggplot2", "igraph", "Hmisc", "ggraph", "RColorBrewer", "scales",
  "circlize", "vegan", "data.table", "mixOmics", "ComplexHeatmap",
  "pheatmap", "ggrepel"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) == 0) {
  message("[OK] All required packages are installed.")
} else {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}
RSCRIPT

if conda run -n "$ENV_NAME" Rscript "$verify_script"; then
  rm -f "$verify_script"
  echo "----------------------------------------------------------------"
  echo "Setup successful."
  echo "Run: conda activate $ENV_NAME"
  echo "----------------------------------------------------------------"
else
  status=$?
  rm -f "$verify_script"
  echo "----------------------------------------------------------------" >&2
  echo "Verification failed. See errors above." >&2
  echo "----------------------------------------------------------------" >&2
  exit "$status"
fi
