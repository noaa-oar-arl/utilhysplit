#!/opt/Tools/anaconda3/bin/python
#mpl.use('Agg')
from optparse import OptionParser
from hcontrol import HycsControl


#pylint:  disable=C0103

parser = OptionParser()

parser.add_option("-i", type="string", dest="fname", default='CONTROL',
                  help="Path and name of input CONTROL file")

(options, args) = parser.parse_args()


control_file = HycsControl(fname=options.fname)
control_file.read()
control_file.rename('CONTROL_annotated')
control_file.write(verbose=True, annotate=True)


