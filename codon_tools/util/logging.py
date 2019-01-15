import logging


# add two custom logging levels
logging.DETAIL = 15
logging.addLevelName(logging.DETAIL, "DETAIL")


def _detail(self, message, *args, **kws):
    if self.isEnabledFor(logging.DETAIL):
        # Yes, logger takes its '*args' as 'args'.
        self._log(logging.DETAIL, message, args, **kws)


logging.OUTPUT = 25
logging.addLevelName(logging.OUTPUT, "OUTPUT")


def _output(self, message, *args, **kws):
    if self.isEnabledFor(logging.OUTPUT):
        # Yes, logger takes its '*args' as 'args'.
        self._log(logging.OUTPUT, message, args, **kws)


logging.Logger.detail = _detail
logging.Logger.output = _output
log_levels = [logging.OUTPUT, logging.INFO, logging.DETAIL, logging.DEBUG]
