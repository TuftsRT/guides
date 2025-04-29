& "$PSScriptRoot/clean.ps1"
if (-not $?) {exit $LASTEXITCODE}
& "$PSScriptRoot/build.ps1"
exit $LASTEXITCODE
