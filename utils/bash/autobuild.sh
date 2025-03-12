#!/bin/bash
root=$(git rev-parse --show-toplevel) || exit 1
sphinx-autobuild --nitpicky --ignore "$root/source/tags" -- "$root/source" "$root/build"
