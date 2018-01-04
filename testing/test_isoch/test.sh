##https://coverage.readthedocs.io/en/coverage-4.1/#quickstart
##"coverage.py is a tool for measuring code coverage of Python programs. It monitors your program noting
##which parts of the code have been exectued, then analyzes the source to identify code that could have been 
##executud but was not."

## Runs the isoch.py program to produce time of arrival maps


rm *png

MDL=$HOME/gitpython/mhysplit

/opt/Tools/anaconda3/bin/python $MDL/isoch.py -icdump.bin

#$MDL/isoch.py -h
#coverage run  $MDL/isoch.py -icdump.bin
#coverage report -m

#coverage run  $MDL/isoch.py -icdump.bin -b30:70:-80:-40 --title 'My test title' --loc 'My test location' -c'CONTROL' 
#coverage report -m


#./pyconcplot -h
#./pyconcplot -icdump_conc

