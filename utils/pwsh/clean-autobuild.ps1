& "$PSScriptRoot/clean.ps1"
if (-not $?) {exit $LASTEXITCODE}
& "$PSScriptRoot/autobuild.ps1"
exit $LASTEXITCODE
