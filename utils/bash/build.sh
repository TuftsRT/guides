#!/bin/bash
root=$(git rev-parse --show-toplevel) || exit $?
sphinx-build --nitpicky "$root/source" "$root/build"
