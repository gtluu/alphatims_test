import argparse
import json
import os


# Parse arguments for CLI usage
def get_args():
    """
    Parse command line parameters, including required, optional, and system parameters.

    :return: Arguments with default or user specified values.
    :rtype: dict
    """
    with open(os.path.join(os.path.realpath(os.path.dirname(__file__)),
                           'parameter_descriptions.json'), 'r') as json_file:
        arg_descriptions = json.loads(json_file.read())

    # Initialize parser.
    parser = argparse.ArgumentParser()

    # Require Arguments
    required = parser.add_argument_group('Required Parameters')
    required.add_argument('--input',
                          help=arg_descriptions['input'],
                          required=True,
                          type=str,
                          nargs='+')

    # Optional Arguments
    optional = parser.add_argument_group('Optional Parameters')
    optional.add_argument('--outdir',
                          help=arg_descriptions['outdir'],
                          default='',
                          type=str)
    optional.add_argument('--mode',
                          help=arg_descriptions['mode'],
                          default='centroid',
                          type=str,
                          choices=['raw', 'centroid', 'profile'])
    optional.add_argument('--compression',
                          help=arg_descriptions['compression'],
                          default='zlib',
                          type=str,
                          choices=['zlib', 'none'])
    optional.add_argument('--ms2_only',
                          help=arg_descriptions['ms2_only'],
                          action='store_true')
    optional.add_argument('--use_raw_calibration',
                          help=arg_descriptions['use_raw_calibration'],
                          action='store_true')
    optional.add_argument('--pressure_compensation_strategy',
                          help=arg_descriptions['pressure_compensation_strategy'],
                          default='global',
                          type=str,
                          choices=['none', 'global', 'frame'])
    optional.add_argument('--exclude_mobility',
                          help=arg_descriptions['exclude_mobility'],
                          action='store_true')
    optional.add_argument('--mz_encoding',
                          help=arg_descriptions['mz_encoding'],
                          default=64,
                          type=int,
                          choices=[32, 64])
    optional.add_argument('--intensity_encoding',
                          help=arg_descriptions['intensity_encoding'],
                          default=64,
                          type=int,
                          choices=[32, 64])
    optional.add_argument('--mobility_encoding',
                          help=arg_descriptions['mobility_encoding'],
                          default=64,
                          type=int,
                          choices=[32, 64])
    optional.add_argument('--barebones_metadata',
                          help=arg_descriptions['barebones_metadata'],
                          action='store_true')
    optional.add_argument('--profile_bins',
                          help=arg_descriptions['profile_bins'],
                          default=0,
                          type=int)
    optional.add_argument('--imzml_mode',
                          help=arg_descriptions['imzml_mode'],
                          default='processed',
                          type=str,
                          choices=['processed', 'continuous'])
    # add another argument for how to handle number of scans per file:
    # Can be minimum, maximum, or mean number of scans per line-scan or a user defined number of scans per line.
    # TODO: Add arg_descriptions for these parameters
    optional.add_argument('--scans_per_line',
                            help=arg_descriptions['scans_per_line'],
                            default='mean',
                            type=str,
                            choices=['minimum', 'maximum', 'mean', 'user_defined'])
    optional.add_argument('--scans_per_line_value', 
                            help=arg_descriptions['scans_per_line_value'],
                            default=0,
                            type=int)

    # System Arguments
    system = parser.add_argument_group('System Parameters')
    system.add_argument('--verbose',
                        help=arg_descriptions['verbose'],
                        action='store_true')

    # Return parser
    arguments = parser.parse_args()
    return vars(arguments)


# Checks to ensure arguments are valid.
def args_check(args):
    """
    Check relevant arguments to ensure user input values are valid.

    :param args: Arguments obtained from timsconvert.arguments.get_args().
    :type args: dict
    """
    # Check if output directory exists and create it if it does not.
    if not os.path.isdir(args['outdir']) and args['outdir'] != '':
        os.makedirs(args['outdir'])
