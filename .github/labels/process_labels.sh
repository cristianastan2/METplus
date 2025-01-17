#!/bin/sh

# Must specify 2 arguments:
#   - username
#   - authorization key

if [[ $# -ne 2 ]]; then
  echo "ERROR: `basename $0` ... must specify the GitHub username and authorization key."
  exit 1
else
  USER=$1
  AUTH=$2
fi

# Script directory
SCRIPT_DIR=`dirname $0`

# List of METplus repositories
REPO_LIST="metplus met metplotpy metcalcpy metdatadb metviewer \
           metexpress metplus-training";

# Build commands to add/update common labels
for REPO in ${REPO_LIST}; do
  echo $REPO
  ${SCRIPT_DIR}/post_patch_labels.sh $USER $AUTH $REPO ${SCRIPT_DIR}/common_labels.txt
done

# Build commands to delete extra labels
for REPO in ${REPO_LIST}; do
  echo $REPO;
  ${SCRIPT_DIR}/delete_labels.sh $USER $AUTH $REPO
done
