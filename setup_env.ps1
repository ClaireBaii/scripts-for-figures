param(
  [string]$EnvFile = "environment.yaml",
  [string]$EnvName = ""
)

$ErrorActionPreference = "Stop"

if ($PSVersionTable.PSVersion.Major -lt 5) {
  throw "This script requires PowerShell 5.0 or later."
}

if (-not (Test-Path $EnvFile)) {
  throw "Environment file not found: $EnvFile"
}

if ([string]::IsNullOrWhiteSpace($EnvName)) {
  $nameLine = Select-String -Path $EnvFile -Pattern '^\s*name\s*:\s*(.+)\s*$' | Select-Object -First 1
  if (-not $nameLine) {
    throw "Failed to parse environment name from $EnvFile"
  }
  $EnvName = $nameLine.Matches[0].Groups[1].Value.Trim()
}

function Resolve-CondaExe {
  if ($env:CONDA_EXE -and (Test-Path $env:CONDA_EXE)) {
    return $env:CONDA_EXE
  }

  $condaCmd = Get-Command conda -ErrorAction SilentlyContinue
  if ($condaCmd) {
    return $condaCmd.Source
  }

  $candidates = @(
    "$env:USERPROFILE\\miniforge3\\condabin\\conda.bat",
    "$env:USERPROFILE\\miniconda3\\condabin\\conda.bat",
    "C:\\ProgramData\\miniforge3\\condabin\\conda.bat",
    "C:\\ProgramData\\miniconda3\\condabin\\conda.bat"
  )

  foreach ($candidate in $candidates) {
    if (Test-Path $candidate) {
      return $candidate
    }
  }

  throw "Conda executable not found. Install Miniforge/Miniconda and run 'conda init powershell'."
}

$CondaExe = Resolve-CondaExe
Write-Host "Using conda: $CondaExe" -ForegroundColor Gray

Write-Host "Updating conda environment '$EnvName' from '$EnvFile'..." -ForegroundColor Cyan
& $CondaExe env update -f $EnvFile --prune
if ($LASTEXITCODE -ne 0) {
  exit $LASTEXITCODE
}

$verifyPath = Join-Path $env:TEMP "verify_packages_${EnvName}_$(Get-Random).R"
$verifyScript = @"
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
"@
Set-Content -Path $verifyPath -Value $verifyScript -Encoding UTF8

try {
  Write-Host "----------------------------------------------------------------" -ForegroundColor Green
  Write-Host "Verifying required R packages in '$EnvName'..." -ForegroundColor Green
  Write-Host "----------------------------------------------------------------" -ForegroundColor Green

  & $CondaExe run -n $EnvName Rscript $verifyPath
  if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
  }

  Write-Host "----------------------------------------------------------------" -ForegroundColor Cyan
  Write-Host "Setup successful." -ForegroundColor Cyan
  Write-Host "Run: conda activate $EnvName" -ForegroundColor White
  Write-Host "----------------------------------------------------------------" -ForegroundColor Cyan
}
finally {
  Remove-Item $verifyPath -ErrorAction SilentlyContinue
}
