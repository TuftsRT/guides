$root = & git rev-parse --show-toplevel
if (-not $?) { exit 1 }
$ProgressPreference = "SilentlyContinue"
$ErrorActionPreference = "SilentlyContinue"
Remove-Item -Path (Join-Path -Path $root -ChildPath "build") -Recurse -Force
Remove-Item -Path (Join-Path -Path $root -ChildPath "jupyter_build") -Recurse -Force
Remove-Item -Path (Join-Path -Path $root -ChildPath "source/tags") -Recurse -Force
