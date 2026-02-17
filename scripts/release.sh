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

time_tag=$(date +%H%M%S)

# commit with date-only message
git commit -m "$date_tag $time_tag"


# push to remote
git push -f

# create release using gh cli
gh release create "$date_tag-$time_tag" --title "$date_tag-$time_tag" --notes "release $date_tag $time_tag"

echo "created release $date_tag $time_tag"

