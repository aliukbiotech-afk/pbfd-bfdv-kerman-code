# Sync local PBFD results to Google Drive with rclone (remote "gdrive:").
# Usage:  powershell -File analysis\sync_to_drive.ps1 [-Source analysis\results] [-Dest PBFD_Kerman_results]
param(
    [string]$Source = (Join-Path $PSScriptRoot "results"),
    [string]$Dest = "PBFD_Kerman_results"
)
$rclone = (Get-Command rclone -ErrorAction SilentlyContinue).Source
if (-not $rclone) {
    $rclone = Get-ChildItem "$env:LOCALAPPDATA\Microsoft\WinGet\Packages" -Recurse -Filter rclone.exe -ErrorAction SilentlyContinue |
              Select-Object -First 1 -ExpandProperty FullName
}
if (-not $rclone) { throw "rclone not found. Install with: winget install Rclone.Rclone" }
if (-not (Test-Path $Source)) { throw "Source folder not found: $Source" }
# 'copy' (not 'sync'): never deletes anything on Drive.
& $rclone copy $Source "gdrive:$Dest" --progress --transfers 4 --exclude "*.ckp.gz" --exclude "__pycache__/**"
if ($LASTEXITCODE -eq 0) { Write-Host "Synced '$Source' -> gdrive:$Dest" } else { Write-Host "rclone exited with code $LASTEXITCODE" }
