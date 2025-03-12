$root = & git rev-parse --show-toplevel
if (-not $?) { exit 1 }
$options = @{
    FilePath = (Get-Command "sphinx-build" -ErrorAction "Stop")
    ArgumentList =
        "--nitpicky",
        (Join-Path -Path $root -ChildPath "source"),
        (Join-Path -Path $root -ChildPath "build")
    NoNewWindow = $true
    Wait = $true
}
Start-Process @options
