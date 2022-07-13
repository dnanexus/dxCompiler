import logging


class Logger(object):

    @classmethod
    def make(cls, verbosity):
        level = {"error": 40, "warning": 30, "info": 20, "debug": 10}[verbosity]
        log_format = '%(asctime)s | %(levelname)8s | %(message)s'
        logger = logging.getLogger(__name__)
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(LogFormatter(log_format))
        logger.addHandler(ch)
        return logger


class LogFormatter(logging.Formatter):
    """
    Logging colored formatter, adapted from https://stackoverflow.com/a/56944256/3638629
    """

    grey = '\x1b[38;21m'
    blue = '\x1b[38;5;39m'
    yellow = '\x1b[38;5;226m'
    red = '\x1b[38;5;196m'
    bold_red = '\x1b[31;1m'
    reset = '\x1b[0m'

    def __init__(self, fmt):
        super().__init__()
        self._fmt = fmt
        self.FORMATS = {
            logging.DEBUG: self.grey + self._fmt + self.reset,
            logging.INFO: self.blue + self._fmt + self.reset,
            logging.WARNING: self.yellow + self._fmt + self.reset,
            logging.ERROR: self.red + self._fmt + self.reset,
            logging.CRITICAL: self.bold_red + self._fmt + self.reset
        }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
