import os
import logging
from multiprocessing import Pool, cpu_count
from timsconvert.constants import VERSION
from timsconvert.data_input import dot_d_detection
from timsconvert.timestamp import get_iso8601_timestamp
from timsconvert_nanodesi.arguments import get_args, args_check
from timsconvert_nanodesi.convert import convert_raw_file, clean_up_logfiles, get_frame_id_for_each_coordinate


def main():
    # Parse arguments.
    args = get_args()
    # Args check.
    args_check(args)
    # Check arguments.
    args['version'] = VERSION

    # Load in input data.
    input_files = []
    for dirpath in args['input']:
        if not dirpath.endswith('.d'):
            input_files = input_files + list(filter(None, dot_d_detection(dirpath)))
        elif dirpath.endswith('.d'):
            if os.path.isdir(dirpath):
                input_files = input_files +[dirpath]
            else:
                logging.info(get_iso8601_timestamp() + ':' + f'{dirpath} does not exist...')
                logging.info(get_iso8601_timestamp() + ':' + 'Skipping...')

    # Get number of scans per line for interpolation.
    frame_ids_at_each_coord = get_frame_id_for_each_coordinate(input_files, args)

    # Convert each sample
    with Pool(processes=cpu_count() - 1) as pool:
        pool_map_input = [(args, infile, line_number+1, frame_ids_at_each_coord[line_number])
                          for line_number, infile in enumerate(input_files)]
        list_of_logfiles = pool.map(convert_raw_file, pool_map_input)
    list_of_logfiles = list(filter(None, list_of_logfiles))

    # Shutdown logger.
    logging.shutdown()

    # Clean up temporary log files.
    clean_up_logfiles(args, list_of_logfiles)


if __name__ == '__main__':
    # Run.
    main()
