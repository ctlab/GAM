#!/bin/bash

DESC_TEMPLATE="GAM/DESCRIPTION.template"
DESC_FILE="GAM/DESCRIPTION"

NEWREV=`git describe --tags --dirty`
VERSION=`git describe --candidates=0 --tags 2>/dev/null`
if [ -z "$VERSION" ] 
then
    VERSION=`git describe --abbrev=0 --tags`-1
fi
DATE=`date '+%Y-%m-%d'`


t=`mktemp`
cat "$DESC_TEMPLATE" |\
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
