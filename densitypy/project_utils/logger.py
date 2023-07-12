import logging
import sys
from datetime import datetime
class ColoredLogger(logging.Formatter):
    """Custom logger formatter to add colors to log messages"""

    COLOR_CODES = {
        "RESET": "\033[0m",
        "INFO": "\033[94m",   # Blue
        "WARNING": "\033[93m",   # Yellow
        "ERROR": "\033[91m",   # Red
    }

    def format(self, record):
        log_time = datetime.utcnow().isoformat() + "Z"
        log_level = record.levelname
        log_message = record.msg

        if log_level == "INFO":
            log_level = self.COLOR_CODES["INFO"] + log_level + self.COLOR_CODES["RESET"]
        elif log_level == "WARNING":
            log_level = self.COLOR_CODES["WARNING"] + log_level + self.COLOR_CODES["RESET"]
        elif log_level == "ERROR":
            log_level = self.COLOR_CODES["ERROR"] + log_level + self.COLOR_CODES["RESET"]

        return f"{log_time} [{log_level}] {log_message}"

# Create logger instance
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create console handler and set formatter
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(ColoredLogger())

# Add console handler to the logger
logger.addHandler(console_handler)
