#!/bin/bash
# Set work directory with argument. 
if [ $# -ne 1 ]; then
  echo "Please specify the working directory." 1>&2
  echo "EXAMPLE:" 1>&2
  echo "    ./cvt2smi.sh souyakuchan_library" 1>&2
  exit 1
fi
WORK_DIR=$1

# Convert sdf to smi.
# It takes about an hour on my PC.
echo 'Convert sdf to smi...'
echo 'Enamine_Advanced_collection'
obabel -isdf $WORK_DIR/Enamine_Advanced_collection.sdf -osmi -O $WORK_DIR/Enamine_Advanced_collection.smi --append idnumber
echo 'Enamine_HTS_collection'
obabel -isdf $WORK_DIR/Enamine_HTS_collection.sdf -osmi -O $WORK_DIR/Enamine_HTS_collection.smi --append idnumber
echo 'Enamine_Premium_collection'
obabel -isdf $WORK_DIR/Enamine_Premium_collection.sdf -osmi -O $WORK_DIR/Enamine_Premium_collection.smi --append idnumber
echo 'UOS_HTS'
obabel -isdf $WORK_DIR/UOS_HTS.sdf -osmi -O $WORK_DIR/UOS_HTS.smi --append ID

echo 'Done.'

