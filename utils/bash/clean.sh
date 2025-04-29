#!/bin/bash
root=$(git rev-parse --show-toplevel) || exit $?
rm -rf "$root/build"
rm -rf "$root/jupyter_execute"
rm -rf "$root/source/tags"
