from timsconvert.parse import parse_lcms_baf, parse_lcms_tsf, parse_lcms_tdf
from timsconvert.classes import TimsconvertBafData, TimsconvertTsfData, TimsconvertTdfData
from timsconvert.timestamp import get_iso8601_timestamp
import os
import sys
import logging
from pyimzml.ImzMLWriter import ImzMLWriter
from pyimzml.compression import NoCompression, ZlibCompression
from pyTDFSDK.util import get_encoding_dtype
# TODO: Emerson: Need to load all input files in order to get all retention times first.
# This lets us get a mask of where to assign spectra for each scan, which aligns spectra spatially.

def write_nanodesi_chunk_to_imzml(data, imzml_file, frame_start, frame_stop, mode, exclude_mobility, profile_bins,
                                  mz_encoding, intensity_encoding, mobility_encoding, line_number, frame_id_for_each_coord):
    """
    Parse and write out a group of spectra to an imzML file from a nano-DESI timsTOF fleX MSI dataset using pyimzML.

    :param data: Object containing raw data information from TDF or TSF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData
    :param imzml_file: Instance of pyimzml.ImzMLWriter.ImzMLWriter for output file.
    :type imzml_file: pyimzml.ImzMLWriter.ImzMLWriter
    :param frame_start: Beginning frame number.
    :type frame_start: int
    :param frame_stop: Ending frame number (non-inclusive).
    :type frame_stop: int
    :param mode: Mode command line parameter, either "profile", "centroid", or "raw".
    :type mode: str
    :param exclude_mobility: Whether to include mobility data in the output files, defaults to None.
    :type exclude_mobility: bool
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param mz_encoding: m/z encoding command line parameter, either "64" or "32".
    :type mz_encoding: int
    :param intensity_encoding: Intensity encoding command line parameter, either "64" or "32".
    :type intensity_encoding: int
    :param mobility_encoding: Mobility encoding command line parameter, either "64" or "32".
    :type mobility_encoding: int
    :param line_number: Line number of the input file.
    :type line_number: int
    :param frame_id_for_each_coord: list containing frame IDs for each coordinate.
    :type frame_id_for_each_coord: list
    """
    # Parse and write TSF data.
    if isinstance(data, TimsconvertTsfData):
        parent_scans, product_scans = parse_lcms_tsf(data,
                                                     frame_start,
                                                     frame_stop,
                                                     mode,
                                                     False,
                                                     profile_bins,
                                                     mz_encoding,
                                                     intensity_encoding)
        # TODO: Test this. It uses a list containing the frame ID to be used in each coordinate.
        for i, scans in parent_scans:
            frame_id = i+frame_start
            # Get index values where frame_id matches the frame_id_for_each_coord list.
            frame_id_matches = [idx for idx, x in enumerate(frame_id_for_each_coord) if x == frame_id]
            for idx in frame_id_matches:
                coord = (line_number, idx)
                imzml_file.addSpectrum(scans.mz_array,
                                       scans.intensity_array,
                                       coord)
    # Parse and write TDF data.
    elif isinstance(data, TimsconvertTdfData):
        parent_scans, product_scans = parse_lcms_tdf(data,
                                                     frame_start,
                                                     frame_stop,
                                                     mode,
                                                     False,
                                                     exclude_mobility,
                                                     profile_bins,
                                                     mz_encoding,
                                                     intensity_encoding,
                                                     mobility_encoding)
        if mode == 'profile':
            exclude_mobility = True
        if not exclude_mobility:
            # TODO: Test this. It uses a list containing the frame ID to be used in each coordinate.
            for i, scans in parent_scans:
                frame_id = i+frame_start
                # Get index values where frame_id matches the frame_id_for_each_coord list.
                frame_id_matches = [idx for idx, x in enumerate(frame_id_for_each_coord) if x == frame_id]
                for idx in frame_id_matches:
                    coord = (line_number, idx)
                    imzml_file.addSpectrum(scans.mz_array,
                                        scans.intensity_array,
                                        coord,
                                        mobilities=scans.mobility_array)
        elif exclude_mobility:
            # TODO: Test this. It uses a list containing the frame ID to be used in each coordinate.
            for i, scans in parent_scans:
                frame_id = i+frame_start
                # Get index values where frame_id matches the frame_id_for_each_coord list.
                frame_id_matches = [idx for idx, x in enumerate(frame_id_for_each_coord) if x == frame_id]
                for idx in frame_id_matches:
                    coord = (line_number, idx)
                    imzml_file.addSpectrum(scans.mz_array,
                                        scans.intensity_array,
                                        coord)
    # Parse and write BAF data.
    elif isinstance(data, TimsconvertBafData):
        parent_scans, product_scans = parse_lcms_baf(data,
                                                     frame_start,
                                                     frame_stop,
                                                     mode,
                                                     False,
                                                     profile_bins,
                                                     mz_encoding,
                                                     intensity_encoding)
        # TODO: Test this. It uses a list containing the frame ID to be used in each coordinate.
        for i, scans in parent_scans:
            frame_id = i+frame_start
            # Get index values where frame_id matches the frame_id_for_each_coord list.
            frame_id_matches = [idx for idx, x in enumerate(frame_id_for_each_coord) if x == frame_id]
            for idx in frame_id_matches:
                coord = (line_number, idx)
                imzml_file.addSpectrum(scans.mz_array,
                                    scans.intensity_array,
                                    coord,
                                    mobilities=scans.mobility_array)



def write_nanodesi_imzml(data, outdir, outfile, mode, exclude_mobility, profile_bins, imzml_mode, mz_encoding,
                         intensity_encoding, mobility_encoding, compression, line_number, frame_id_for_each_coord, chunk_size=10):
    """
    Parse and write out spectra to an imzML file from a nano-DESI timsTOF fleX MSI dataset using pyimzML.

    :param data: Object containing raw data information from TDF or TSF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData
    :param outdir: Output directory path that was specified from the command line parameters or the original input
        file path if no output directory was specified.
    :type outdir: str
    :param outfile: The original input filename if no output filename was specified.
    :type outfile: str
    :param mode: Mode command line parameter, either "profile", "centroid", or "raw".
    :type mode: str
    :param exclude_mobility: Whether to include mobility data in the output files, defaults to None.
    :type exclude_mobility: bool
    :param profile_bins: Number of bins to bin spectrum to.
    :type profile_bins: int
    :param imzml_mode: Whether to export spectra in "processed" (individual m/z and intensity arrays per pixel) or
        "continuous" mode (single m/z array for the entire dataset, individual intensity arrays per pixel).
    :type imzml_mode: str
    :param mz_encoding: m/z encoding command line parameter, either "64" or "32".
    :type mz_encoding: int
    :param intensity_encoding: Intensity encoding command line parameter, either "64" or "32".
    :type intensity_encoding: int
    :param mobility_encoding: Mobility encoding command line parameter, either "64" or "32".
    :type mobility_encoding: int
    :param compression: Compression command line parameter, either "zlib" or "none".
    :type compression: str
    :param chunk_size: Number of MS1 spectra that to be used when subsetting dataset into smaller groups to pass onto
        timsconvert.write.write_lcms_chunk_to_mzml() for memory efficiency; larger chunk_size requires more memory
        during conversion.
    :type chunk_size: int
    """
    # Set polarity for run in imzML.
    polarity = list(set(data.analysis['Frames']['Polarity'].values.tolist()))
    if len(polarity) == 1 and polarity[0] == '+':
        polarity = 'positive'
    elif len(polarity) == 1 and polarity[0] == '-':
        polarity = 'negative'
    else:
        polarity = None

    if data.analysis['GlobalMetadata']['SchemaType'] == 'TSF' and mode == 'raw':
        logging.info(get_iso8601_timestamp() + ':' + 'TSF file detected. Only export in profile or centroid mode are '
                                                     'supported. Defaulting to centroid mode.')

    # Get compression type object.
    if compression == 'zlib':
        compression_object = ZlibCompression()
    elif compression == 'none':
        compression_object = NoCompression()

    if data.analysis['GlobalMetadata']['SchemaType'] == 'TSF':
        writer = ImzMLWriter(os.path.join(outdir, outfile),
                             polarity=polarity,
                             mode=imzml_mode,
                             spec_type=mode,
                             mz_dtype=get_encoding_dtype(mz_encoding),
                             intensity_dtype=get_encoding_dtype(intensity_encoding),
                             mz_compression=compression_object,
                             intensity_compression=compression_object,
                             include_mobility=False)
    elif data.analysis['GlobalMetadata']['SchemaType'] == 'TDF':
        if mode == 'profile':
            exclude_mobility = True
            logging.info(
                get_iso8601_timestamp() + ':' + 'Export of ion mobility data is not supported for profile mode data...')
            logging.info(get_iso8601_timestamp() + ':' + 'Exporting without ion mobility data...')
        if not exclude_mobility:
            writer = ImzMLWriter(os.path.join(outdir, outfile),
                                 polarity=polarity,
                                 mode=imzml_mode,
                                 spec_type=mode,
                                 mz_dtype=get_encoding_dtype(mz_encoding),
                                 intensity_dtype=get_encoding_dtype(intensity_encoding),
                                 mobility_dtype=get_encoding_dtype(mobility_encoding),
                                 mz_compression=compression_object,
                                 intensity_compression=compression_object,
                                 mobility_compression=compression_object,
                                 include_mobility=True)
        elif exclude_mobility:
            writer = ImzMLWriter(os.path.join(outdir, outfile),
                                 polarity=polarity,
                                 mode=imzml_mode,
                                 spec_type=mode,
                                 mz_dtype=get_encoding_dtype(mz_encoding),
                                 intensity_dtype=get_encoding_dtype(intensity_encoding),
                                 mz_compression=compression_object,
                                 intensity_compression=compression_object,
                                 include_mobility=False)

    with writer as imzml_file:
        chunk = 0
        frames = data.analysis['Frames']['Id'].to_list()
        while chunk + chunk_size + 1 <= len(frames):
            chunk_list = []
            for i, j in zip(frames[chunk:chunk + chunk_size], frames[chunk + 1: chunk + chunk_size + 1]):
                chunk_list.append((int(i), int(j)))
            logging.info(get_iso8601_timestamp() +
                         ':' +
                         'Parsing and writing Frame ' +
                         str(chunk_list[0][0]) +
                         ' from ' +
                         data.analysis['GlobalMetadata']['SampleName'] +
                         '...')
            for frame_start, frame_stop in chunk_list:
                write_nanodesi_chunk_to_imzml(data,
                                              imzml_file,
                                              frame_start,
                                              frame_stop,
                                              mode,
                                              exclude_mobility,
                                              profile_bins,
                                              mz_encoding,
                                              intensity_encoding,
                                              mobility_encoding,
                                              line_number, 
                                              frame_id_for_each_coord)
                sys.stdout.write(get_iso8601_timestamp() +
                                 ':' +
                                 data.source_file.replace('/', '\\') +
                                 ':Progress:' +
                                 str(round((frame_start / data.analysis['Frames'].shape[0]) * 100)) +
                                 '%\n')
            chunk += chunk_size
        else:
            chunk_list = []
            for i, j in zip(frames[chunk:-1], frames[chunk + 1:]):
                chunk_list.append((int(i), int(j)))
            chunk_list.append((j, data.analysis['Frames'].shape[0] + 1))
            logging.info(get_iso8601_timestamp() +
                         ':' +
                         'Parsing and writing Frame ' +
                         str(chunk_list[0][0]) +
                         ' from ' +
                         data.analysis['GlobalMetadata']['SampleName'] +
                         '...')
            for frame_start, frame_stop in chunk_list:
                write_nanodesi_chunk_to_imzml(data,
                                              imzml_file,
                                              frame_start,
                                              frame_stop,
                                              mode,
                                              exclude_mobility,
                                              profile_bins,
                                              mz_encoding,
                                              intensity_encoding,
                                              mobility_encoding,
                                              line_number, 
                                              frame_id_for_each_coord)
                sys.stdout.write(get_iso8601_timestamp() +
                                 ':' +
                                 data.source_file.replace('/', '\\') +
                                 ':Progress:100%\n')
    logging.info(
        get_iso8601_timestamp() + ':' + 'Finished writing to .imzML file ' + os.path.join(outdir, outfile) + '...')
