import sys
import argparse
from fragment_analyzer import Fsa
import logging
import yaml
from yaml.loader import SafeLoader
from fragment_analyzer.method.method import create_method_obj
from pathlib import Path
from fragment_analyzer.plot.plot import plot_trace

def create_logger(logpath: Path | None = None, loglevel: str = "INFO") -> logging.Logger:
    """

    :param loglevel: string, default 'INFO', possible values: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    :param logpath: string or None, default None
    :return: logging.Logger object
    """

    level = logging.INFO

    match loglevel:
        case "DEBUG":
            level = logging.DEBUG
        case "INFO":
            level = logging.INFO
        case "WARNING":
            level = logging.WARNING
        case "ERROR":
            level = logging.ERROR
        case "CRITICAL":
            level = logging.CRITICAL

    FORMAT = "%(asctime)s;%(levelname)s;%(message)s"

    handlers = [logging.StreamHandler()]
    if logpath:
        handlers.append(logging.FileHandler(logpath, 'w+'))

    logging.basicConfig(handlers=handlers, level=level, format=FORMAT)
    return logging.getLogger('fraggler_log')


def get_arguments(preset_methods: list) -> argparse.Namespace:
    """
    
    :param preset_methods: list of built-in methods apply
    :return: argparse.Namespace object
    """

    preset_loglevels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

    parser = argparse.ArgumentParser(
        prog='fraggler',
        description='Automatic fragment analyzer')

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--infile', type=Path, help='input fragment analysis file (.fsa)')
    input_group.add_argument('-d', '--indir', type=Path, help='dir containing fragment analysis files (dir/*.fsa)')

    method_group = parser.add_mutually_exclusive_group(required=True)
    method_group.add_argument('-m',
                              '--method',
                              type=str,
                              help='select a method from available options',
                              choices=preset_methods)

    method_group.add_argument('-mp', '--method-path', type=Path, help='path to method yaml file (.yaml)')

    parser.add_argument('-l', '--log-file-path',  type=Path, help='save log to logfile', required=False)
    parser.add_argument('-x', '--log-level', default='INFO',
                        help='set log level (default: INFO)', choices=preset_loglevels, required=False)

    args = parser.parse_args()
    print(args)

    if args.infile is not None:
        if not args.infile.is_file():
            parser.print_help()
            exit(1)


    return parser.parse_args()

def main():

    methods_dirpath = Path("fragment_analyzer/config/methods")
    method_filepaths = methods_dirpath.glob("*.yaml")
    methods_dict = {}
    for m in method_filepaths:
        methods_dict[m.name.split('.')[0]] = m

    args = get_arguments(list(methods_dict.keys()))

    method_path = args.method_path if args.method_path else methods_dict[args.method]
    logger = create_logger(args.log_file_path, args.log_level)
    logger.info(f'Method path {method_path} selected')

    try:
        with open(method_path) as f:
            method = yaml.load(f, Loader=SafeLoader)
            logger.info(f'Method path {method_path} loaded')
    except Exception as e:
        logger.error(f'Method path {method_path} cannot be loaded')
        for row in str(e).split('\n'):
            logger.error(f"Stack trace:     {row}")

        exit(1)

    try:
        logger.info(f"Name of selected method in yaml is {method['name']}")
    except Exception as e:
        logger.error(f"Method name in from method file cannot be printed. The selected file is probably not a method file")
        for row in str(e).split('\n'):
            logger.error(f"Stack trace:     {row}")

        exit(1)

    method_obj = create_method_obj(method, logger=logger)
    logger.info(f'Method object for {method_obj.name} created ')

    fsa_obj = Fsa(file=args.infile, method=method_obj, correct_baseline=True)
    logger.info(f'FSA object for {fsa_obj.filename} using method {method_obj.name} created ')

    plot_trace(fsa_obj, "ladder", raw=True)


if __name__ == '__main__':
    main()

