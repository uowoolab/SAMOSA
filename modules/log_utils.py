from __future__ import annotations

# handle differing python versions (< 3.8, & > 3.9) of CSD python API
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

import logging

from enum import Enum
from multiprocessing import Queue

from logging.handlers import QueueHandler, QueueListener


# logging level for type hinting
class LogLevel(Enum):
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL


LogLevelType = Literal[
    LogLevel.DEBUG, LogLevel.INFO, LogLevel.WARNING, LogLevel.ERROR, LogLevel.CRITICAL
]


# main process logging
def main_logging_setup(log_queue: Queue, log_level: LogLevelType) -> logging.Logger:
    """
    Set up logging in the main process.

    Parameters:
        log_queue (multiprocessing.Queue): shared queue for child processes logging.
        log_level (logging.LOGLEVEL): desired level of logging information.

    Returns:
        logger (logging.Logger): main process logging.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # remove duplicate handlers
    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    s_handler = logging.StreamHandler()
    s_handler.setFormatter(get_formatter())
    root_logger.addHandler(s_handler)

    f_handler = logging.FileHandler(filename="samosa.log", mode="a")
    f_handler.setFormatter(get_formatter())
    root_logger.addHandler(f_handler)

    # Configure a queue listener
    listener = QueueListener(log_queue, s_handler)
    listener.start()
    return root_logger, listener


#
def child_logging_setup(
    log_queue: Queue, process_name: str, log_level: LogLevelType
) -> logging.Logger:
    """
    Set up logging in the main process.

    Parameters:
        log_queue (multiprocessing.Queue): shared queue for child processes logging.
        process (str): child process name.
        log_level (logging.LOGLEVEL): desired level of logging information.

    Returns:
        logger (logging.Logger): main process logging.
    """
    logger = logging.getLogger(process_name)
    logger.propagate = False
    logger.setLevel(log_level)

    # Add a QueueHandler to send logs to the main process
    queue_handler = QueueHandler(log_queue)
    if not any(isinstance(h, QueueHandler) for h in logger.handlers):
        logger.addHandler(queue_handler)
    if not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
        file_handler = logging.FileHandler(filename="samosa.log", mode="a")
        file_handler.setFormatter(get_formatter())
        logger.addHandler(file_handler)
    return logger


# log message format
def get_formatter() -> logging.Formatter:
    return logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s (%(funcName)s) | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# get logger
def get_logger(name) -> logging.Logger:
    return logging.getLogger(name)

