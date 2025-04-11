#!/bin/bash
root=$(git rev-parse --show-toplevel) || exit $?
sphinx-autobuild --nitpicky --ignore "$root/source/tags" -- "$root/source" "$root/build"
