import os
import datetime
import matplotlib

class ColorMaker:
    def __init__(self, cname, nvals, ctype="hex", transparency="C8"):
        """
        cname : name of matplotlib colormap
        nvals : number of color values
        ctype : if 'hex' returns 8 digit hexidecimal with
                transparancy.
        transparency : str: transparency value to use in hexidecimal.
        """
        self.transparency = transparency
        self.clist = []
        self.ctype = ctype
        self.get_cmap(cname, nvals)

    def __call__(self):
        """
        Returns:
        list of nvals colors equally spaced throughout the colormap.
        and in hexidecimal format.
        """
        return self.clist

    def rgb_to_hex(self, rgb):
        """
        convert from rgb to hexidecimal.
        """

        def subfunc(val):
            rval = hex(int(val * 255)).replace("0x", "").upper()
            if len(rval) == 1:
                rval = "0" + rval
            return rval

        hval = [subfunc(x) for x in list(rgb)]
        if self.transparency:
            return "{}{}{}{}".format(self.transparency, hval[0], hval[1], hval[2])
        else:
            return "{}{}{}".format(hval[0], hval[1], hval[2])

    def get_cmap(self, cname="viridis", nvals=10):
        cmap = matplotlib.cm.get_cmap(cname)
        cvals = cmap.N
        cspace = int(cvals / nvals)
        if self.ctype == "hex":
            self.clist = [self.rgb_to_hex(cmap(x)) for x in range(0, cvals, cspace)]
        else:
            self.clist = [cmap(x) for x in range(0, cvals, cspace)]

        # for iii in range(0,cvals,cspace):
        #    if ctype == 'hex':
        #        self.clist.append(self.rgb_to_hex(cmap(iii)))
        #    else:
        #        self.clist.append(cmap[iii])


