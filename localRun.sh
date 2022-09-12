#!/bin/bash
python3 ./tools/DatafedVerify.py
retVal=$?
if [ $retVal -eq 0 ];
then
cd runs
for f in *.sh
do
bash "$f"
done
cd ..
echo "Runs completed"
else
echo "Datafed is not configured properly"
fi
