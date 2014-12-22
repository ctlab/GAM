#!/bin/bash

DESC_FILE="GAM/DESCRIPTION"

VERSION=`git describe --dirty --tags`
DATE=`date '+%Y-%m-%d'`


t=`mktemp`
cat "$DESC_FILE" |\
    sed "s/^Revision:.*$/Revision: $NEWREV/" |\
    sed "s/^Date:.*$/Date: $DATE/" |\
    sed "s/^Version:.*$/Version: $VERSION/" |\
    cat - > "$t"


if ! cmp --quiet "$t" "$DESC_FILE"
then
    echo "Updating description"
    cp "$t" "$DESC_FILE"
fi
rm "$t"
