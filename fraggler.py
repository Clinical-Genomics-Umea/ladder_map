import argparse
from fragment_analyzer import Fsa
import logging
import yaml
from yaml.loader import SafeLoader
# from fragment_analyzer.validators.method_validator import validator
from fragment_analyzer.validators.method_validation import validate_method
from pathlib import Path


def get_arguments(methods: list):
    # methods_path = Path("fragment_analyzer/config/methods")
    # method_objects = methods_path.glob("*.yaml")
    # method_dict = {}
    # for m in method_objects:
    #     method_dict[m.name.split('.')[0]] = m

    parser = argparse.ArgumentParser(
        prog='fraggler',
        description='Automatic fragment analyzer')

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--infile', type=str, help='Input fragment analysis file (fsa)')
    input_group.add_argument('-d', '--indir', type=str, help='Dir containing fragment analysis files (fsa)')

    method_group = parser.add_mutually_exclusive_group(required=True)
    method_group.add_argument('-m',
                              '--method',
                              type=str,
                              help='Input fragment analysis file (fsa)',
                              choices=methods)

    method_group.add_argument('-mp', '--method-path', type=str, help='Path to method yaml file')

    parser.add_argument('-l', '--log-file-path',  type=str, help='Save log to logfile')

    return parser.parse_args()

def main():

    #
    # parser = argparse.ArgumentParser(
    #     prog='fraggler',
    #     description='Automatic fragment analyzer')
    #
    # input_group = parser.add_mutually_exclusive_group(required=True)
    # input_group.add_argument('-i', '--infile', type=str, help='Input fragment analysis file (fsa)')
    # input_group.add_argument('-d', '--indir', type=str, help='Dir containing fragment analysis files (fsa)')
    #
    # method_group = parser.add_mutually_exclusive_group(required=True)
    # method_group.add_argument('-m',
    #                           '--method',
    #                           type=str,
    #                           help='Input fragment analysis file (fsa)',
    #                           choices=methods_dict.keys())
    #
    # method_group.add_argument('-mp', '--method-path', type=str, help='Path to method yaml file')
    #
    # parser.add_argument('-l', '--log-file-path',  type=str, help='Save log to logfile')

    methods_path = Path("fragment_analyzer/config/methods")
    method_objects = methods_path.glob("*.yaml")
    methods_dict = {}
    for m in method_objects:
        methods_dict[m.name.split('.')[0]] = m

    args = get_arguments(list(methods_dict.keys()))
    method_path = args.method_path if args.method_path else methods_dict[args.method]

    with open(method_path) as f:
        method = yaml.load(f, Loader=SafeLoader)
        logging.info(f'method path {method_path} loaded')

    method_model = validate_method(method)

    print(method_model)

    # if validator(method):
    #     logging.info(f'method definition ok')

    # fsa_obj = Fsa(file=args.infile, method=method, correct_baseline=True)

    # print(repr(fsa_obj))



if __name__ == '__main__':
    main()

