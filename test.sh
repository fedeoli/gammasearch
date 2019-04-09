echo "TEST: clear_env.py"
python clear_env.py
echo "PASSED"

echo "TEST: recall_inaf.py"
python recall_inaf.py 220 45 4 10
echo "PASSED"

echo "TEST: zipper.py"
python zipper.py Generated_Maps Generated_Events 220 45
echo "PASSED"

echo "TEST: python_launch"
python -W ignore python_launch.py zipdir/total/Map_1.fits single.log b 2.3 0.8 3 23 3.9 66 y
echo "PASSED"

echo "TEST: data_analysis"
python -W ignore data_analysis.py zipdir/Map.log zipdir/total final.log b 2.3 0.8 3 23 3.9 66
echo "PASSED"

echo "TEST: general_analysis"
python -W ignore general_analysis.py zipdir/total final.log b 2.3 0.8 3 23 3.9 66
echo "PASSED"

echo "CLEARING THE FILESYSTEM"
python clear_env.py
