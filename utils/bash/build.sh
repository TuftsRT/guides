#!/bin/bash
root=$(git rev-parse --show-toplevel) || exit 1
sphinx-build --nitpicky "$root/source" "$root/build"
