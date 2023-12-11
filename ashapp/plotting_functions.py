#!/opt/Tools/anaconda3/envs/hysplit/bin/python
# -----------------------------------------------------------------------------
# Air Resources Laboratory
#
# plotting_functions.py
#
# -----------------------------------------------------------------------------
#
#
# -----------------------------------------------------------------------------

# from abc mport ABC, abstractmethod
import datetime
import glob
import logging
import os
import subprocess
import sys
import zipfile

from utilvolc.runhelper import Helper
from utilhysplit.plotutils.concplotutil import ConcplotColors

logger = logging.getLogger(__name__)


def cleanup(self):
    # TO DO test and modify
    psfiles = glob.glob("*.ps")
    for psf in psfiles:
        Helper.remove(psf)
    return -1


def check_output_type(outputname):
    temp = outputname.split(".")
    if "ps" in temp[-1].lower():
        return "+g0"
    elif "html" in temp[-1].lower():
        return "+g1"
    elif "svg" in temp[-1].lower():
        return "+g1"
    else:
        return "+g0"


def create_parxplot(inp, pardump_filename, outputname):
    # 11/13/2020 do not create kml for pardump.
    # files too large and not useful.
    maparg = "-j" + os.path.join(inp["MAP_DIR"], "arlmap")
    JOBID = inp["jobid"]
    convert_exe = inp["CONVERT_EXE"]
    ghostscript_exe = inp["GHOSTSCRIPT_EXE"]
    resolution = inp["graphicsResolution"]

    output_type = check_output_type(outputname)
    # gisoption = self.inp['gisOption']
    gisoption = 0
    c = [
        os.path.join(inp["HYSPLIT_DIR"], "exec", "parxplot"),
        "-i" + pardump_filename,
        "-k1",
        output_type,
        "-o" + outputname,
        # "-p{}".format(JOBID),
        "-a{}".format(gisoption),
        maparg,
    ]
    logger.warning(" ".join(c))
    Helper.execute(c)
    rval = check_and_convert(
        outputname,
        output_type,
        JOBID,
        convert_exe,
        ghostscript_exe,
        final_type="gif",
        resolution=resolution,
    )
    return rval


def check_and_convert(
    outputname,
    intial_type,
    JOBID,
    convert_exe,
    ghostscript_exe,
    final_type="gif",
    resolution=80,
):
    if not os.path.exists(outputname):
        logger.debug(outputname)
        logger.error(
            "******************************************************************************"
        )
        logger.error(
            "The model was not able to create {} plot for job {}.".format(
                outputname, JOBID
            )
        )
        logger.error(
            "******************************************************************************"
        )
        rval = (False, JOBID, "FAILED")
    else:
        if intial_type == "+g0":
            convert_ps_to_image(
                convert_exe,
                ghostscript_exe,
                outputname,
                outputname.replace("ps", final_type),
                resolution=resolution,
            )
        rval = (True, JOBID, "PASSED")
    return rval


def create_massloading_plot(inp, inputname, outputname, stage=0, conc_multiplier=1):
    # create dispersion png plot using python concplot.
    # Units in g/m2
    gisopt = inp["gisOption"]
    mapopt = inp["mapBackground"]
    mapdir = inp["MAP_DIR"]
    hdir = inp["HYSPLIT_DIR"]
    JOBID = inp["jobid"]
    convert_exe = inp["CONVERT_EXE"]
    ghostscript_exe = inp["GHOSTSCRIPT_EXE"]
    resolution = inp["graphicsResolution"]

    clrs = ConcplotColors()
    # want in grams. output in mg.
    conc_multiplier = conc_multiplier / 1000.0
    contours = "+".join(
        [
            "50:High:{}".format(clrs.get("red")),
            "10:High:{}".format(clrs.get("magenta")),
            "5.0:Medium:{}".format(clrs.get("orange")),
            "2.0:Medium:{}".format(clrs.get("yellow")),
            "1.0:Low:{}".format(clrs.get("green")),
            "0.1::{}".format(clrs.get("tan")),
            "0.01::{}".format(clrs.get("grey")),
        ]
    )
    logger.info("Creating column mass loading graphics for job {}.".format(JOBID))
    if mapopt == "terrain":
        maparg = "--street-map=0"
    elif mapopt == "toner":
        maparg = "--street-map=1"
    else:
        maparg = "-j" + os.path.join(mapdir, "arlmap")

    fortran_concplot = True
    output_type = check_output_type(outputname)
    if fortran_concplot:
        # c = [self.inp['PYTHON_EXE'],
        #     os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),
        c = [
            os.path.join(hdir, "exec", "concplot"),
            "-i" + inputname,
            "-c4",
            output_type,  # svg or ps
            "-v" + contours,
            "-e4",  # mass loading
            "-d4",  # over all levels
            "-ug",  # units of m2 rather than m3.
            #'+n',
            # this will change otput name of kml file.
            # "-p" + str(JOBID),
            "-x{:2.2e}".format(conc_multiplier),
            "-o" + outputname,
            "-s0",  # sum species
            "".join(["-j", hdir, "/graphics/arlmap"]),
            #'-m{}'.format(self.inp['mapProjection']),
            #'-z{}'.format(self.inp['zoomFactor']),
            #'-l1',
            maparg,
            "-a{}".format(gisopt),
            #'-g0:{}'.format(self.inp['spatialPlotRadius']),
            "-k1",
        ]
    # else
    #    c = [self.inp['PYTHON_EXE'],
    #         os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),

    logger.info("Executing concplot " + " ".join(c))
    Helper.execute(c)
    # outputname = outputname.replace("ps", JOBID)
    rval = check_and_convert(
        outputname,
        output_type,
        JOBID,
        convert_exe,
        ghostscript_exe,
        final_type="png",
        resolution=resolution,
    )
    # TODO - put this in a different function?
    #        test that it works.
    # make the kmz file.
    # if the -p option used need to use the get_gelabel_filename method.
    # gelist = self.make_gelabel(self.filelocator.get_gelabel_filename)
    # self.make_gelabel(self.filelocator.get_basic_gelabel_filename)
    # fnames = [self.filelocator.get_basic_kml_filename()]
    # fnames = ["HYSPLIT_ps.kml"]
    # fnames = [self.filelocator.get_basic_kml_filename()]
    # logger.info("KML FILE NAME {}".format(fnames[0]))
    # fnames.extend(gelist)
    # logger.info("GELIST FILE NAME {}".format(fnames[1]))
    # self.generate_kmz(
    #    fnames, self.filelocator.get_kmz_filename(stage="massload")
    # )
    # remove the kml and gelabel files
    # for fn in fnames:
    #    Helper.remove(fn)
    return rval


def make_gelabel(hysplitdir, filenamer, JOBID):
    # TO DO - make it work.

    # the get_basic_gelabel_filename method should be passed if the -p
    # option was not used and files have form GELABEL_??_ps.txt.
    # otherwise the get_gelabel_filename method should be passed.

    # run gelabel executable.
    # this assumes there is a gelabel.txt file created by concplot.
    # gelabel executable creates ps files of the legend from the txt file with

    c = [
        os.path.join(hysplitdir, "exec", "gelabel"),
    ]
    # check to see if -p option is needed.
    if filenamer(frame=1, ptype="txt").split(".")[0][-2:] != "ps":
        c.append("-p{}".format(JOBID))
    logger.debug(" ".join(c))
    Helper.execute(c)

    iii = 1
    # the ps files are then converted to gif files.
    gelabel = filenamer(frame=iii, ptype="ps")
    logger.debug(gelabel + " NAME ")
    gelist = []
    # while os.path.isfile(gelabel):
    #    gelabelgif = filenamer(frame=iii, ptype="gif")
    #    gelist.append(gelabelgif)
    #    logger.debug(gelabel + " NAME " + gelabelgif)
    #    self.convert_ps_to_image(gelabel, gelabelgif, resolution=80)
    #    iii += 1
    #    gelabel = filenamer(frame=iii, ptype="ps")
    # if iii > 100:
    #    break
    return gelist


def create_concentration_plot(inp, inputname, outputname, conc_multiplier=1):
    """
    currently generates kml file which overwrites the massloading kml file.
    but does not generate a kmz file.

    took out the -p option which makes filename.jobid
    """
    JOBID = inp["jobid"]
    hdir = inp["HYSPLIT_DIR"]
    mapopt = inp["mapBackground"]
    mapdir = inp["MAP_DIR"]
    convert_exe = inp["CONVERT_EXE"]
    ghostscript_exe = inp["GHOSTSCRIPT_EXE"]
    output_type = check_output_type(outputname)
    resolution = inp["graphicsResolution"]

    # create dispersion png plot using python concplot.
    clrs = ConcplotColors()
    contours = "+".join(
        [
            "1000:High:{}".format(clrs.get("purple")),
            "100:High:{}".format(clrs.get("red")),
            "10:High:{}".format(clrs.get("magenta")),
            "5.0:Medium:{}".format(clrs.get("orange")),
            "2.0:Medium:{}".format(clrs.get("yellow")),
            "0.2:Low:{}".format(clrs.get("green")),
            "0.02:None:{}".format(clrs.get("tan")),
            "0.002:None:{}".format(clrs.get("grey")),
        ]
    )
    logger.info("Creating concentration graphics for job {}.".format(JOBID))
    gisopt = 3  # self.inp['gisOption']

    if mapopt == "terrain":
        maparg = "--street-map=0"
    elif mapopt == "toner":
        maparg = "--street-map=1"
    else:
        maparg = "-j" + os.path.join(mapdir, "arlmap")

    fortran_concplot = True
    if fortran_concplot:
        # c = [self.inp['PYTHON_EXE'],
        #     os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),
        c = [
            os.path.join(hdir, "exec", "concplot"),
            "-i" + inputname,
            "-c4",
            output_type,  # output in ps or svg
            "-v" + contours,
            "-umg",  # output in mg
            #'+n',
            # "-p" + str(JOBID),
            "-x{:2.2e}".format(conc_multiplier),
            "-o" + outputname,
            "-s0",  # sum species
            "".join(["-j", hdir, "/graphics/arlmap"]),
            #'-m{}'.format(self.inp['mapProjection']),
            #'-z{}'.format(self.inp['zoomFactor']),
            #'-l1',
            # maparg,
            "-a{}".format(gisopt),
            #'-g0:{}'.format(self.inp['spatialPlotRadius']),
            "-k1",
        ]
    # else
    #    c = [self.inp['PYTHON_EXE'],
    #         os.path.join(self.inp['HYSPLIT_DIR'], 'exec', 'concplot.py'),
    logger.debug("Executing concplot " + " ".join(c))
    Helper.execute(c)
    rval = check_and_convert(
        outputname,
        output_type,
        JOBID,
        convert_exe,
        ghostscript_exe,
        final_type="png",
        resolution=resolution,
    )
    return rval


def create_montage_page(inp, inputnamelist, outputname, stage, iii, jjj):
    # TODO modify and test.
    """
    flin : function which generates filename
    flout : function which generates filename
    jjj : frame number of output
    iii : frame number of first input file
    stage : input for flin and flout
    """
    done = False
    # self.filelocator.get_concentration_montage_filename(stage,frame=jjj)
    # jlist, levlist = set_qva_levels()
    # rlist = []
    montage = inp["CONVERT_EXE"].replace("convert", "montage")
    c = [montage]
    for inputname in inputnamelist:
        lev = 0
        # for lev in levlist:
        #    self.filelocator.get_concplot_filename(stage=stage,frame=iii)
        if not os.path.exists(inputname):
            done = True
            logger.debug("file does not exist {}".format(inputname))
            break
        c.append("-label '{}'".format(lev))
        c.append("-pointsize 20")
        c.append("-trim {}".format(inputname))
        iii += 1
    if done:
        return None, done, iii
    c.extend(["-geometry 200x200", "-tile 2x6", outputname])
    # logger.debug("Create montage: " + " ".join(c))
    # Sonny - not sure why I need shell=True?
    # subprocess.Popen(' '.join(c), shell=True)
    Helper.execute_with_shell(c)
    return outputname, done, iii


def create_concentration_montage(self, stage):
    flout = self.filelocator.get_concentration_montage_filename
    flin = self.filelocator.get_concplot_filename
    file_addlist = [self.filelocator.get_massloading_filename]
    file_addlist.append(self.filelocator.get_parxplot_filename)
    stagein = [stage, stage, stage]
    self.create_montage_pdf(stagein, flin, flout, file_addlist)


def create_montage_pdf(self, stage, flin, flout, file_addlist):
    """
    stage : list of stages.
    flin : function which gennerates file names
    flout : function which gennerates file names
    file_addlist: list of function which generate file names
    """
    # TO DO - a lot of loops and very clunky.
    # Assumes 4 levels. To change need to change 'tile'
    logger.debug("Creating pdf montage")
    done = False
    montage_list = []
    iii = 0
    jjj = 0
    while not done:
        outputname, done, iii = self.create_montage_page(
            flin, flout, stage[0], iii, jjj
        )
        # add the other files to the list for the pdf file.
        nnn = 1
        for fl in file_addlist:
            newname = fl(stage=stage[nnn], frame=jjj, ptype="gif")
            if os.path.isfile(newname):
                montage_list.append(newname)
            else:
                logger.debug("file not found {}".format(newname))
            nnn += 1

        if outputname:
            montage_list.append(outputname)
        if jjj >= 100:
            done = True
        jjj += 1
    if montage_list:
        c = [self.inp["CONVERT_EXE"]]
        c.extend(montage_list)
        c.append(self.filelocator.get_totalpdf_filename(stage=stage[0]))
        logger.debug("MONTAGE PDF {}".format(" ".join(c)))
        Helper.execute_with_shell(c)


def convert_ps_to_image(
    convert_exe, ghostscript_exe, ps_filename, gif_filename, resolution=80, trim=True
):
    logger.debug("Creating images from ps {}".format(ps_filename))
    # if resolution == 0:
    #    resolution = inp["graphicsResolution"]

    if not os.path.exists(ps_filename):
        logger.warn(
            "Postscript file {} does not exist. Image file {} will not be created".format(
                ps_filename, gif_filename
            )
        )
        return
    if trim:
        c = [
            convert_exe,
            "+adjoin",
            "-",
            "-trim",
            "+repage",
            "GIF:{}".format(gif_filename),
        ]
    else:
        c = [
            convert_exe,
            "+adjoin",
            "-",
            "+repage",
            "GIF:{}".format(gif_filename),
        ]
    p1 = subprocess.Popen(
        [
            ghostscript_exe,
            "-r{}".format(resolution),
            "-dTextAlphaBits=4",
            "-dGraphicsAlphaBits=4",
            "-dNOPAUSE",
            "-dSAFER",
            "-sDEVICE=pnmraw",
            "-q",
            "-sOutputFile=-",
            ps_filename,
            "-c",
            "quit",
        ],
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
    )
    p2 = subprocess.Popen(
        c,
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
    )
    p1.stdout.close()  # allow p1 to receive a SIGPIPE if p2 exists.
    stdoutdata, stderrdata = p2.communicate()


def convert_ps_to_pdf(inp, ps_filename, pdf_filename):
    if not os.path.exists(ps_filename):
        logger.warn(
            "Postscript file {} does not exist. PDF file {} will not be created".format(
                ps_filename, pdf_filename
            )
        )
        return

    cproc = [
        inp["GHOSTSCRIPT_EXE"],
        "-q",
        "-dBATCH",
        "-dSAFER",
        "-dMaxBitmap=500000000",
        "-dNOPAUSE",
        "-dAlignToPixels=0",
        "-dEPSCrop",
        "-sDEVICE=pdfwrite",
        "-sOutputFile={}".format(pdf_filename),
        "-f{}".format(ps_filename),
    ]
    logger.info("Creating PDF " + " ".join(cproc))
    Helper.execute(cproc)


def get_maptext_info(inp, metfilefinder):
    maptexthash = {}
    rstr = "HYSPLIT dispersion run. "
    maptexthash["run_description"] = rstr

    metfiles = metfilefinder.find(inp["start_date"], inp["durationOfSimulation"])
    if metfiles:
        metfiles = list(zip(*metfiles))[1]
    else:
        metfiles = "unknown"

    maptexthash["infoc"] = ",".join(metfiles)

    return maptexthash


def create_maptext(inp, maptexthash, conc_multiplier, ash_reduction):
    # TO DO - write meteorology
    # vertical distribution.
    # MER and m63.
    JOBID = inp["jobid"]
    descripline = maptexthash["run_description"]

    now = datetime.datetime.utcnow().strftime("%m/%d/%Y %H:%M UTC")

    emission_rate = conc_multiplier()  # gives rate in mg/hour
    emission_rate = "{:2.1e} kg/h".format(emission_rate / 1e6)

    name = inp["VolcanoName"]

    emission = float(inp["emissionHours"])
    ehours = int(emission)
    eminutes = round((emission - int(emission)) * 60)

    m63 = ash_reduction

    lat = inp["latitude"]
    lon = inp["longitude"]
    alt = inp["bottom"]
    top = inp["top"]
    sp8 = "        "
    estart = inp["start_date"].strftime("%Y %m %d %H:%M UTC")

    # metid = self.inp["meteorologicalData"]
    # start_lon
    # start_lat
    # data[0] #alert type
    # vaac
    # hgts[0]
    # data[7]  #Alert detection time
    jobidline = "Job ID: {}   Job Name: {}   Job Completion: {}\n".format(
        JOBID, inp["jobname"], now
    )
    volcanoline = "Volcano: {}{} lat: {} lon: {}{}Hgt: {} to {} m\n".format(
        name, sp8, lat, lon, sp8, alt, top
    )
    poll = "Pollutant: Ash \n"
    ## TO DO calculate release quantity
    release_a = "Start: {}  Duration: {} h {} min\n".format(estart, ehours, eminutes)
    release_b = "Release: {}    m63: {:1g}\n".format(emission_rate, m63)
    info_a = "Vertical distribution: {}{}  GSD: {}{}  Particles: {}\n".format(
        "uniform", sp8, "default", sp8, "20,000"
    )
    info_b = "Meteorology: {}\n".format(inp["meteorologicalData"])
    info_c = maptexthash["infoc"]
    owner = "Run by: {}\n".format(inp["owner"])

    with open(inp["WORK_DIR"] + "MAPTEXT.CFG", "w") as fid:
        fid.write("  \n")
        fid.write("  \n")
        fid.write(jobidline)
        fid.write("  \n")
        fid.write(descripline)
        fid.write("  \n")
        fid.write("  \n")
        fid.write(volcanoline)
        fid.write(poll)
        fid.write(release_a)
        fid.write(release_b)
        fid.write(info_a)
        fid.write(info_b)
        fid.write(info_c)
        fid.write(owner)
        # current = dtime.utcnow()
        # fid.write('Job Start: '+current.strftime('%B %d, %Y')+' at '+current.strftime('%H:%M:%S')+' UTC \n')
        fid.write("\n")
        fid.close()


# def generate_shape_file(self):
# no shape files generated now.


def generate_kmz(hysplit_dir, kml_filenames, kmz_filename, compresslevel):
    for f in kml_filenames:
        if not os.path.exists(f):
            logger.warn(
                "KML file {} does not exist. Google Earth file {} will not be created".format(
                    f, kmz_filename
                )
            )
            return

    files = [
        os.path.join(hysplit_dir, "guicode", "noaa_google.gif"),
        os.path.join(hysplit_dir, "guicode", "logocon.gif"),
        os.path.join(hysplit_dir, "graphics", "blueball.png"),
        os.path.join(hysplit_dir, "graphics", "greenball.png"),
        os.path.join(hysplit_dir, "graphics", "redball.png"),
    ]
    files += kml_filenames

    with zipfile.ZipFile(
        kmz_filename, "w", compresslevel=compresslevel
    ) as z:
        for f in files:
            if os.path.exists(f):
                bn = os.path.basename(f)
                z.write(f, arcname=bn)


def create_zipped_up_file(inp, filename, files):
    if not files:
        return False
    logger.debug("{} files to be zipped {}".format(filename, "\n".join(files)))
    with zipfile.ZipFile(
        filename, "w", compresslevel=compresslevel
    ) as z:
        for f in files:
            if os.path.exists(f):
                z.write(f)
    return True
