
import logging

PARANOID = 5
name = "PARANOID"
class logger_with_extra_verbose_mode(logging.Logger):
    def __init__(self,name):
        logging.Logger.__init__(self,name)

    def paranoid(self,text,*args):
        self.log(PARANOID,text,*args)

if not logging.getLevelName(PARANOID) == name:
    logging.addLevelName(name,PARANOID)
logging.setLoggerClass(logger_with_extra_verbose_mode)

