make -C ~/hemelb/Code
rm -r results

./run.sh

diff -q CleanImages results/Images

#if [ ! -d results/Snapshots/converted ]; then
# mkdir results/Snapshots/converted
#fi

#ls -1 results/Snapshots/ | grep ".dat" | while read FILE
#do
# ./XdrToText results/Snapshots/$FILE results/Snapshots/converted/$FILE
#done

./NumericalComparison CleanSnapshots results/Snapshots
