from optparse import OptionParser
from models/hcontrol import HycsControl


#pylint:  disable=C0103

parser = OptionParser()

parser.add_option("-i", type="string", dest="fname", default='./CONTROL',
                  help="Path and name of input CONTROL file")
parser.add_option("-h", action='store_true', dest="help", default=False,
                  help="print out help message")



(options, args) = parser.parse_args()
if options.help:
    hstr = '-i path/name of CONTROL file to check. Default is ./CONTROL \n'
    hstr +='Uses the HycsControl class to read in the file.'
    hstr +='outputs an annotated file called CONTROL_annotated'

else:
    control_file = HycsControl(fname=options.fname)
    control_file.read()
    control_file.rename('CONTROL_annotated')
    control_file.write(verbose=True, annotate=True)


