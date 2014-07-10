#!/bin/bash

DESC_FILE="GAM/DESCRIPTION"

NEWREV=`hg id -i`
NEWREV_NUM=`hg id -n | sed "s/+/.5/"`
VERSION=`hg id -t -r 'ancestors(.) and tag()'`
TAGREV=`hg id -i -r 'ancestors(.) and tag()'`
if [ "$TAGREV" != "$NEWREV" ]
then
    VERSION="${VERSION}-${NEWREV_NUM}"
fi
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
