$root = & git rev-parse --show-toplevel
if (-not $?) { exit 1 }
$options = @{
    FilePath = (Get-Command "sphinx-autobuild" -ErrorAction "Stop")
    ArgumentList =
        "--nitpicky",
        "--ignore $(Join-Path -Path $root -ChildPath "source/tags")",
        "--",
        (Join-Path -Path $root -ChildPath "source"),
        (Join-Path -Path $root -ChildPath "build")
    NoNewWindow = $true
    Wait = $true
}
Start-Process @options
