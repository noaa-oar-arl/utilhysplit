
class ConcplotColors:
    def __init__(self):
        colorhash = {}
        colorhash["yellow"] = "242236051"
        colorhash["orange"] = "235137052"
        colorhash["red"] = "161024014"
        colorhash["blue"] = "070051242"
        colorhash["green"] = "147219121"
        colorhash["magenta"] = "194056143"
        colorhash["purple"] = "107023156"
        colorhash["cyan"] = "075201199"
        colorhash["grey"] = "150150150"
        colorhash["tan"] = "163145131"
        self.colorhash = colorhash

    def get(self, color):
        return self.colorhash[color]

