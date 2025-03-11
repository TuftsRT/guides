& "$PSScriptRoot/clean.ps1"
if (-not $?) {exit 1}
& "$PSScriptRoot/autobuild.ps1"
#exit $?
