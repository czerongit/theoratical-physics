#!/bin/bash
# commit staged files and create github release

set -e

# check if any files are staged
staged_count=$(git diff --cached --name-only | wc -l | tr -d ' ')

if [ "$staged_count" -eq 0 ]; then
    echo "no staged files, exiting"
    exit 1
fi

# generate date string MMDDYYYY
date_tag=$(date +%m%d%Y)

# commit with date-only message
git commit -m "$date_tag"

# push to remote
git push

# create release using gh cli
gh release create "$date_tag" --title "$date_tag" --notes "release $date_tag"

echo "created release $date_tag"

