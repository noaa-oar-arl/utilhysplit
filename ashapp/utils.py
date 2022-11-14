import logging
import sys

MESSAGE_FORMAT_INFO = "%(asctime)s.%(msecs)03d %(levelname)s - %(message)s"
MESSAGE_FORMAT_DEBUG = \
    "%(asctime)s.%(msecs)03d %(levelname)s %(name)s - %(message)s"
    
def setup_logger(level=logging.INFO):
    logging.basicConfig(
            stream=sys.stdout,
            level=level,
            format=MESSAGE_FORMAT_INFO,
            datefmt="%H:%M:%S")
