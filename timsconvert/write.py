from timsconvert.parse import (parse_lcms_baf, parse_lcms_tsf, parse_lcms_tdf, parse_maldi_tsf, parse_maldi_tdf,
                               parse_maldi_tdf_iprm, parse_maldi_plate_map)
from timsconvert.classes import TimsconvertBafData, TimsconvertTsfData, TimsconvertTdfData
from timsconvert.constants import VERSION, MSMS_TYPE_CATEGORY, INSTRUMENT_SOURCE_TYPE
from timsconvert.timestamp import get_iso8601_timestamp
import os
import sys
import logging
import copy
import numpy as np
from psims.mzml import MzMLWriter
from pyteomics import mgf
from pyimzml.ImzMLWriter import ImzMLWriter
from pyimzml.compression import NoCompression, ZlibCompression
from pyTDFSDK.util import get_encoding_dtype


def write_mzml_metadata(data, writer, infile, mode, ms2_only, barebones_metadata):
    """
    Write metadata to mzML file using psims. Includes spectral metadata, source files, software list, instrument
    configuration, and data processing. Note that TIMSCONVERT is not included as a data processing step due to a lack
    CV param and compatibility with downstream data analysis software.

    :param data: Object containing raw data information from TDF, TSF, or BAF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData |
        timsconvert.classes.TimsconvertBafData
    :param writer: Instance of psims.mzml.MzMLWriter for output file.
    :type writer: psims.mzml.MzMLWriter
    :param infile: Input file path to be used for source file metadata.
    :type infile: str
    :param mode: Mode command line parameter, either "profile", "centroid", or "raw".
    :type mode: str
    :param ms2_only: Whether to include MS1 data in the output files.
    :type ms2_only: bool
    :param barebones_metadata: If True, omit software and data processing metadata in the resulting mzML files. Used
        for compatibility with downstream analysis software that does not have support for newer CV params or
        UserParams.
    :type barebones_metadata: bool
    """
    # Basic file descriptions.
    file_description = []
    if isinstance(data, TimsconvertBafData):
        metadata_key = 'Properties'
    elif isinstance(data, TimsconvertTsfData) or isinstance(data, TimsconvertTdfData):
        metadata_key = 'GlobalMetadata'

    # Add spectra level and centroid/profile status.
    if isinstance(data, TimsconvertBafData):
        ms_levels = list(set(data.analysis['AcquisitionKeys']['MsLevel'].values.tolist()))
        ms_levels = [int(i) for i in ms_levels]
        if 0 in ms_levels:
            file_description.append('MS1 spectrum')
        if 1 in ms_levels:
            file_description.append('MSn spectrum')
    elif isinstance(data, TimsconvertTsfData) or isinstance(data, TimsconvertTdfData):
        ms_levels = list(set(data.analysis['Frames']['MsMsType'].values.tolist()))
        ms_levels = [int(i) for i in ms_levels]
        ms_levels_tmp = []
        for i in ms_levels:
            if i in MSMS_TYPE_CATEGORY['ms1']:
                ms_levels_tmp.append('MS1 spectrum')
            elif i in MSMS_TYPE_CATEGORY['ms2']:
                ms_levels_tmp.append('MSn spectrum')
        ms_levels_tmp = list(set(ms_levels_tmp))
        for i in ms_levels_tmp:
            file_description.append(i)
    else:
        if not ms2_only:
            file_description.append('MS1 spectrum')
            file_description.append('MSn spectrum')
        elif ms2_only:
            file_description.append('MSn spectrum')
    if mode == 'raw' or mode == 'centroid':
        file_description.append('centroid spectrum')
    elif mode == 'profile':
        file_description.append('profile spectrum')

    # Source file
    sf = writer.SourceFile(os.path.split(infile)[0],
                           os.path.split(infile)[1],
                           id=os.path.splitext(os.path.split(infile)[1])[0])
    writer.file_description(file_contents=file_description, source_files=sf)

    # Add list of software.
    if not barebones_metadata:
        acquisition_software = {'id': 'micrOTOFcontrol',
                                'version': data.analysis[metadata_key]['AcquisitionSoftwareVersion'],
                                'params': ['micrOTOFcontrol', ]}
        tdf_sdk_software = {'id': 'TDF-SDK',
                            'version': '2.21.0.4',
                            'params': ['Bruker software']}
        psims_software = {'id': 'psims-writer',
                          'version': '1.3.4',
                          'params': ['python-psims', ]}
        timsconvert_software = {'id': 'timsconvert',
                                'version': VERSION,
                                'params': ['timsconvert', ]}
        writer.software_list([acquisition_software,
                              tdf_sdk_software,
                              psims_software,
                              timsconvert_software])

    # Instrument configuration.
    # Instrument source, analyzer, and detector are all hard coded to timsTOF hardware and does not allow for non-stock
    # configurations.
    inst_count = 1
    if data.analysis[metadata_key]['InstrumentSourceType'] in INSTRUMENT_SOURCE_TYPE.keys() \
            and 'MaldiApplicationType' not in data.analysis[metadata_key].keys():
        source = writer.Source(inst_count,
                               [INSTRUMENT_SOURCE_TYPE[data.analysis[metadata_key]['InstrumentSourceType']]])
    elif 'MaldiApplicationType' in data.analysis[metadata_key].keys():
        source = writer.Source(inst_count, ['matrix-assisted laser desorption ionization'])
    # Analyzer and detector hard coded for timsTOF fleX
    inst_count += 1
    analyzer = writer.Analyzer(inst_count, ['quadrupole', 'time-of-flight'])
    inst_count += 1
    detector = writer.Detector(inst_count, ['microchannel plate detector', 'photomultiplier'])
    # Get instrument serial number.
    serial_number = data.analysis[metadata_key]['InstrumentSerialNumber']
    # Get instrument name based on GlobalMetadata or Properties table.
    if not barebones_metadata:
        instrument_name = data.analysis[metadata_key]['InstrumentName'].strip().lower()
        if instrument_name == 'timstof':
            instrument_name = 'timsTOF'
        elif instrument_name == 'timstof pro':
            instrument_name = 'timsTOF Pro'
        elif instrument_name == 'timstof pro 2':
            instrument_name = 'timsTOF Pro 2'
        elif instrument_name == 'timstof flex':
            instrument_name = 'timsTOF fleX'
        elif instrument_name == 'timstof scp':
            instrument_name = 'timsTOF SCP'
        elif instrument_name == 'timstof ht':
            instrument_name = 'Bruker Daltonics timsTOF series'  # placeholder since HT doesn't have CV param
        elif instrument_name == 'timstof ultra':
            instrument_name = 'timsTOF Ultra'
        params = [instrument_name, {'instrument serial number': serial_number}]
    else:
        params = [{'instrument serial number': serial_number}]
    inst_config = writer.InstrumentConfiguration(id='instrument',
                                                 component_list=[source, analyzer, detector],
                                                 params=params)
    writer.instrument_configuration_list([inst_config])

    # Data processing element.
    if not barebones_metadata:
        proc_methods = [writer.ProcessingMethod(order=1,
                                                software_reference='psims-writer',
                                                params=['Conversion to mzML']),
                        writer.ProcessingMethod(order=2,
                                                software_reference='TDF-SDK',
                                                params=['Conversion to mzML']),
                        writer.ProcessingMethod(order=3,
                                                software_reference='timsconvert',
                                                params=['Conversion to mzML'])]
        processing = writer.DataProcessing(proc_methods, id='exportation')
        writer.data_processing_list([processing])


def get_spectra_count(data):
    """
    Calculate the predicted number of spectra to be written.

    :param data: Object containing raw data information from TDF, TSF, or BAF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData |
        timsconvert.classes.TimsconvertBafData
    :return: Number of expected spectra.
    :rtype: int
    """
    if isinstance(data, TimsconvertTdfData):
        ms1_count = data.analysis['Frames'][data.analysis['Frames']['MsMsType'] == 0]['MsMsType'].values.size
        if 'Precursors' in data.analysis.keys():
            ms2_count = len(list(filter(None, data.analysis['Precursors']['MonoisotopicMz'].values)))
        # Set ms2_count to 0 if precursors table is not found.
        else:
            ms2_count = 0
        ms_count = ms1_count + ms2_count
    elif isinstance(data, TimsconvertTsfData):
        ms_count = data.analysis['Frames'].shape[0]
    elif isinstance(data, TimsconvertBafData):
        ms1_count = data.analysis['Spectra'][data.analysis['Spectra']['AcquisitionKey'] ==
                                             1]['AcquisitionKey'].values.size
        ms2_count = data.analysis['Spectra'][data.analysis['Spectra']['AcquisitionKey'] ==
                                             2]['AcquisitionKey'].values.size
        ms_count = ms1_count + ms2_count
    return ms_count


def update_spectra_count(outdir, outfile, num_of_spectra, scan_count):
    """
    Calculate the actual number of spectra that were written to the output mzML file. Update is needed to prevent
    counting emtpy spectra that were omitted from the output. Performs this by iterating over tmp mzML file that was
    written and replacing the tmp file with a final mzML file with a true spectra count.

    :param outdir: Output directory path that was specified from the command line parameters or the original input
        file path if no output directory was specified.
    :type outdir: str
    :param outfile: The original input filename if no output filename was specified.
    :type outfile: str
    :param num_of_spectra: Number of spectra that was calculated for the current file being converted using
        timsconvert.write.get_spectra_count().
    :type num_of_spectra: int
    :param scan_count: Final true count for the number of spectra from the current file being converted.
    :type scan_count: int
    """
    if os.path.exists(os.path.join(outdir, outfile)):
        os.remove(os.path.join(outdir, outfile))
    with open(os.path.splitext(os.path.join(outdir, outfile))[0] + '_tmp.mzML', 'r') as in_stream, \
            open(os.path.join(outdir, outfile), 'w') as out_stream:
        for line in in_stream:
            out_stream.write(line.replace(
                '      <spectrumList count="' + str(num_of_spectra) + '" defaultDataProcessingRef="exportation">',
                '      <spectrumList count="' + str(scan_count) + '" defaultDataProcessingRef="exportation">'))
    os.remove(os.path.splitext(os.path.join(outdir, outfile))[0] + '_tmp.mzML')


def write_ms1_spectrum(writer, data, scan, mz_encoding, intensity_encoding, mobility_encoding, compression, title=None):
    """
    Write an MS1 spectrum to an mzML file using psims.

    :param writer: Instance of psims.mzml.MzMLWriter for output file.
    :type writer: psims.mzml.MzMLWriter
    :param data: Object containing raw data information from TDF, TSF, or BAF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData |
        timsconvert.classes.TimsconvertBafData
    :param scan: Object containing spectrum metadata and data arrays.
    :type scan: pyBaf2Sql.classes.BafSpectrum | pyTDFSDK.classes.TsfSpectrum | pyTDFSDK.classes.TdfSpectrum
    :param mz_encoding: m/z encoding command line parameter, either "64" or "32".
    :type mz_encoding: int
    :param intensity_encoding: Intensity encoding command line parameter, either "64" or "32".
    :type intensity_encoding: int
    :param mobility_encoding: Mobility encoding command line parameter, either "64" or "32".
    :type mobility_encoding: int
    :param compression: Compression command line parameter, either "zlib" or "none".
    :type compression: str
    :param title: Spectrum title to be used for MALDI data, defaults to None.
    :type title: str | None
    """
    if isinstance(data, TimsconvertBafData):
        metadata_key = 'Properties'
    elif isinstance(data, TimsconvertTsfData) or isinstance(data, TimsconvertTdfData):
        metadata_key = 'GlobalMetadata'
    # Build params list for spectrum.
    params = [scan.scan_type,
              {'ms level': scan.ms_level},
              {'total ion current': scan.total_ion_current},
              {'base peak m/z': scan.base_peak_mz},
              # {'base peak intensity': scan.base_peak_intensity},
              ({'name': 'base peak intensity',
                'unit_name': 'number of detector counts',
                'value': scan.base_peak_intensity}),
              {'highest observed m/z': scan.high_mz},
              {'lowest observed m/z': scan.low_mz}]
    if 'MaldiApplicationType' in data.analysis[metadata_key].keys():
        params.append({'maldi spot identifier': scan.coord})
        params.append({'spectrum title': title})
    if scan.ms2_no_precursor:
        params.append({'collision energy': scan.collision_energy})

    if scan.mobility_array is not None:
        other_arrays = [({'name': 'mean inverse reduced ion mobility array',
                          'unit_name': 'volt-second per square centimeter'},
                         scan.mobility_array)]
        # other_arrays = [('mean inverse reduced ion mobility array', scan.mobility_array)]
    else:
        other_arrays = None

    encoding_dict = {'m/z array': get_encoding_dtype(mz_encoding),
                     'intensity array': get_encoding_dtype(intensity_encoding)}
    if other_arrays is not None:
        encoding_dict['mean inverse reduced ion mobility array'] = get_encoding_dtype(mobility_encoding)

    writer.write_spectrum(scan.mz_array,
                          scan.intensity_array,
                          id='scan=' + str(scan.scan_number),
                          polarity=scan.polarity,
                          centroided=scan.centroided,
                          scan_start_time=scan.retention_time,
                          other_arrays=other_arrays,
                          params=params,
                          encoding=encoding_dict,
                          compression=compression)


def write_ms2_spectrum(writer, data, scan, mz_encoding, intensity_encoding, mobility_encoding, compression,
                       parent_scan=None, title=None):
    """
    Write an MS/MS spectrum to an mzML file using psims.

    :param writer: Instance of psims.mzml.MzMLWriter for output file.
    :type writer: psims.mzml.MzMLWriter
    :param data: Object containing raw data information from TDF, TSF, or BAF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData |
        timsconvert.classes.TimsconvertBafData
    :param scan: Object containing spectrum metadata and data arrays.
    :type scan: pyBaf2Sql.classes.BafSpectrum | pyTDFSDK.classes.TsfSpectrum | pyTDFSDK.classes.TdfSpectrum
    :param mz_encoding: m/z encoding command line parameter, either "64" or "32".
    :type mz_encoding: int
    :param intensity_encoding: Intensity encoding command line parameter, either "64" or "32".
    :type intensity_encoding: int
    :param mobility_encoding: Mobility encoding command line parameter, either "64" or "32".
    :type mobility_encoding: int
    :param compression: Compression command line parameter, either "zlib" or "none".
    :type compression: str
    :param parent_scan: Object containing standard spectrum data for parent scan used to link MS/MS spectrum with
        an MS1 spectrum.
    :type parent_scan: pyBaf2Sql.classes.BafSpectrum | pyTDFSDK.classes.TsfSpectrum | pyTDFSDK.classes.TdfSpectrum
    :param title: Spectrum title to be used for MALDI data, defaults to None.
    :type title: str | None
    """
    if isinstance(data, TimsconvertBafData):
        metadata_key = 'Properties'
    elif isinstance(data, TimsconvertTsfData) or isinstance(data, TimsconvertTdfData):
        metadata_key = 'GlobalMetadata'
    # Build params list for spectrum.
    params = [scan.scan_type,
              {'ms level': scan.ms_level},
              {'total ion current': scan.total_ion_current}]
    if 'MaldiApplicationType' in data.analysis[metadata_key].keys():
        params.append({'spectrum title': title})
    if hasattr(scan, 'base_peak_mz') and hasattr(scan, 'base_peak_intensity'):
        params.append({'base peak m/z': scan.base_peak_mz})
        # params.append({'base peak intensity': scan.base_peak_intensity})
        params.append(({'name': 'base peak intensity',
                        'unit_name': 'number of detector counts',
                        'value': scan.base_peak_intensity}))
    if hasattr(scan, 'high_mz') and hasattr(scan, 'low_mz'):
        params.append({'highest observed m/z': scan.high_mz})
        params.append({'lowest observed m/z': scan.low_mz})

    if scan.mobility_array is not None:
        other_arrays = [({'name': 'mean inverse reduced ion mobility array',
                          'unit_name': 'volt-second per square centimeter'},
                         scan.mobility_array)]
        # other_arrays = [('mean inverse reduced ion mobility array', scan['mobility_array'])]
    else:
        other_arrays = None

    encoding_dict = {'m/z array': get_encoding_dtype(mz_encoding),
                     'intensity array': get_encoding_dtype(intensity_encoding)}
    if other_arrays is not None:
        encoding_dict['mean inverse reduced ion mobility array'] = get_encoding_dtype(mobility_encoding)

    # Build precursor information dict.
    precursor_info = {'mz': scan.selected_ion_mz,
                      'activation': ['collision-induced dissociation',  # hard coded to either CID or isCID
                                     {'collision energy': scan.collision_energy}],
                      'isolation_window_args': {'target': scan.target_mz,
                                                'upper': scan.isolation_upper_offset,
                                                'lower': scan.isolation_lower_offset},
                      'params': []}

    if scan.selected_ion_intensity is not None:
        precursor_info['intensity'] = scan.selected_ion_intensity
    if scan.selected_ion_mobility is not None:
        precursor_info['params'].append({'inverse reduced ion mobility': scan.selected_ion_mobility})
    if scan.selected_ion_ccs is not None:
        precursor_info['params'].append({'collisional cross sectional area': scan.selected_ion_ccs})
    if scan.charge_state is not None and not np.isnan(scan.charge_state):
        if int(scan.charge_state) != 0:
            precursor_info['charge'] = scan.charge_state

    if parent_scan is not None:
        precursor_info['spectrum_reference'] = 'scan=' + str(parent_scan.scan_number)

    # Write out MS2 spectrum.
    writer.write_spectrum(scan.mz_array,
                          scan.intensity_array,
                          id='scan=' + str(scan.scan_number),
                          polarity=scan.polarity,
                          centroided=scan.centroided,
                          scan_start_time=scan.retention_time,
                          other_arrays=other_arrays,
                          params=params,
                          precursor_information=precursor_info,
                          encoding=encoding_dict,
                          compression=compression)


def write_lcms_chunk_to_mzml(data, writer, frame_start, frame_stop, scan_count, mode, ms2_only, exclude_mobility,
                             profile_bins, mz_encoding, intensity_encoding, mobility_encoding, compression):
    """
    Parse and write out a group of spectra to an mzML file from an LC-MS(/MS) dataset using psims.

    :param data: Object containing raw data information from TDF, TSF, or BAF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData |
        timsconvert.classes.TimsconvertBafData
    :param writer: Instance of psims.mzml.MzMLWriter for output file.
    :type writer: psims.mzml.MzMLWriter
    :param frame_start: Beginning frame number.
    :type frame_start: int
    :param frame_stop: Ending frame number (non-inclusive).
    :type frame_stop: int
    :param scan_count: Current count for the number of spectra from the current file that have been converted.
    :type scan_count: int
    :param mode: Mode command line parameter, either "profile", "centroid", or "raw".
    :type mode: str
    :param ms2_only: Whether to include MS1 data in the output files.
    :type ms2_only: bool
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
    :param compression: Compression command line parameter, either "zlib" or "none".
    :type compression: str
    :return: Updated count for the number of spectra from the current file that have been converted.
    :rtype: int
    """
    # Parse TDF data
    if isinstance(data, TimsconvertTdfData):
        parent_scans, product_scans = parse_lcms_tdf(data,
                                                     frame_start,
                                                     frame_stop,
                                                     mode,
                                                     ms2_only,
                                                     exclude_mobility,
                                                     profile_bins,
                                                     mz_encoding,
                                                     intensity_encoding,
                                                     mobility_encoding)
    # Parse TSF data
    elif isinstance(data, TimsconvertTsfData):
        parent_scans, product_scans = parse_lcms_tsf(data,
                                                     frame_start,
                                                     frame_stop,
                                                     mode,
                                                     ms2_only,
                                                     profile_bins,
                                                     mz_encoding,
                                                     intensity_encoding)
    # Parse BAF data
    elif isinstance(data, TimsconvertBafData):
        parent_scans, product_scans = parse_lcms_baf(data,
                                                     frame_start,
                                                     frame_stop,
                                                     mode,
                                                     ms2_only,
                                                     profile_bins,
                                                     mz_encoding,
                                                     intensity_encoding)

    # Write MS1 parent scans.
    if not ms2_only and product_scans != []:
        for parent in parent_scans:
            products = [i for i in product_scans if i.parent_frame == parent.frame]
            # Set params for scan.
            scan_count += 1
            parent.scan_number = scan_count
            write_ms1_spectrum(writer, data, parent, mz_encoding, intensity_encoding, mobility_encoding, compression)
            # Write MS2 Product Scans
            for product in products:
                scan_count += 1
                product.scan_number = scan_count
                write_ms2_spectrum(writer, data, product, mz_encoding, intensity_encoding, mobility_encoding,
                                   compression, parent_scan=parent)
    elif ms2_only or parent_scans == []:
        for product in product_scans:
            scan_count += 1
            product.scan_number = scan_count
            write_ms2_spectrum(writer, data, product, mz_encoding, intensity_encoding, mobility_encoding, compression)
    elif not product_scans:
        for scan in parent_scans:
            scan_count += 1
            scan.scan_number = scan_count
            if scan.ms_level == 1:
                write_ms1_spectrum(writer, data, scan, mz_encoding, intensity_encoding, mobility_encoding, compression)
            elif scan.ms_level == 2:
                if scan.ms2_no_precursor:
                    write_ms1_spectrum(writer, data, scan, mz_encoding, intensity_encoding, mobility_encoding,
                                       compression)
                else:
                    write_ms2_spectrum(writer, data, scan, mz_encoding, intensity_encoding, mobility_encoding,
                                       compression)
    return scan_count


def write_lcms_mzml(data, infile, outdir, outfile, mode, ms2_only, exclude_mobility, profile_bins, mz_encoding,
                    intensity_encoding, mobility_encoding, compression, barebones_metadata, chunk_size=10, gui=False):
    """
    Parse and write out spectra to an mzML file from an LC-MS(/MS) dataset using psims.

    :param data: Object containing raw data information from TDF, TSF, or BAF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData |
        timsconvert.classes.TimsconvertBafData
    :param infile: Input file path to be used for source file metadata.
    :type infile: str
    :param outdir: Output directory path that was specified from the command line parameters or the original input
        file path if no output directory was specified.
    :type outdir: str
    :param outfile: The original input filename if no output filename was specified.
    :type outfile: str
    :param mode: Mode command line parameter, either "profile", "centroid", or "raw".
    :type mode: str
    :param ms2_only: Whether to include MS1 data in the output files.
    :type ms2_only: bool
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
    :param compression: Compression command line parameter, either "zlib" or "none".
    :type compression: str
    :param barebones_metadata: If True, omit software and data processing metadata in the resulting mzML files. Used
        for compatibility with downstream analysis software that does not have support for newer CV params or
        UserParams.
    :type barebones_metadata: bool
    :param chunk_size: Number of MS1 spectra that to be used when subsetting dataset into smaller groups to pass onto
        timsconvert.write.write_lcms_chunk_to_mzml() for memory efficiency; larger chunk_size requires more memory
        during conversion.
    :type chunk_size: int
    :param gui: Whether function is running from GUI or command line. Toggles sys.stdout messages.
    :type gui: bool
    """
    if isinstance(data, TimsconvertBafData):
        frames_key = 'Spectra'
        metadata_key = 'Properties'
    elif isinstance(data, TimsconvertTsfData) or isinstance(data, TimsconvertTdfData):
        frames_key = 'Frames'
        metadata_key = 'GlobalMetadata'

    # Initialize mzML writer using psims.
    logging.info(get_iso8601_timestamp() + ':' + 'Initializing mzML Writer...')
    writer = MzMLWriter(os.path.splitext(os.path.join(outdir, outfile))[0] + '_tmp.mzML', close=True)

    with writer:
        # Begin mzML with controlled vocabularies (CV).
        logging.info(get_iso8601_timestamp() + ':' + 'Initializing controlled vocabularies...')
        writer.controlled_vocabularies()

        # Start write acquisition, instrument config, processing, etc. to mzML.
        logging.info(get_iso8601_timestamp() + ':' + 'Writing mzML metadata...')
        write_mzml_metadata(data, writer, infile, mode, ms2_only, barebones_metadata)

        logging.info(
            get_iso8601_timestamp() + ':' + 'Writing data to .mzML file ' + os.path.join(outdir, outfile) + '...')
        # Parse chunks of data and write to spectrum elements.
        with writer.run(id='run',
                        instrument_configuration='instrument',
                        start_time=data.analysis[metadata_key]['AcquisitionDateTime']):
            scan_count = 0
            # Count number of spectra in run.
            logging.info(get_iso8601_timestamp() + ':' + 'Calculating number of spectra...')
            num_of_spectra = get_spectra_count(data)
            with writer.spectrum_list(count=num_of_spectra):
                chunk = 0
                # Write data in chunks of chunk_size.
                while chunk + chunk_size + 1 <= len(data.ms1_frames):
                    chunk_list = []
                    for i, j in zip(data.ms1_frames[chunk: chunk + chunk_size],
                                    data.ms1_frames[chunk + 1: chunk + chunk_size + 1]):
                        chunk_list.append((int(i), int(j)))
                    logging.info(get_iso8601_timestamp() +
                                 ':' +
                                 'Parsing and writing Frame ' +
                                 str(chunk_list[0][0]) +
                                 ' from ' +
                                 os.path.split(data.source_file.replace('/', '\\'))[1] +
                                 '...')
                    for frame_start, frame_stop in chunk_list:
                        scan_count = write_lcms_chunk_to_mzml(data,
                                                              writer,
                                                              frame_start,
                                                              frame_stop,
                                                              scan_count,
                                                              mode,
                                                              ms2_only,
                                                              exclude_mobility,
                                                              profile_bins,
                                                              mz_encoding,
                                                              intensity_encoding,
                                                              mobility_encoding,
                                                              compression)
                        if gui:
                            sys.stdout.write(get_iso8601_timestamp() +
                                             ':' +
                                             data.source_file.replace('/', '\\') +
                                             ':Progress:' +
                                             str(round((frame_start / data.analysis[frames_key].shape[0]) * 100)) +
                                             '%\n')
                    chunk += chunk_size
                # Last chunk may be smaller than chunk_size
                else:
                    chunk_list = []
                    for i, j in zip(data.ms1_frames[chunk:-1], data.ms1_frames[chunk + 1:]):
                        chunk_list.append((int(i), int(j)))
                    chunk_list.append((j, data.analysis[frames_key].shape[0] + 1))
                    logging.info(get_iso8601_timestamp() +
                                 ':' +
                                 'Parsing and writing Frame ' +
                                 str(chunk_list[0][0]) +
                                 ' from ' +
                                 os.path.split(data.source_file.replace('/', '\\'))[1] +
                                 '...')
                    for frame_start, frame_stop in chunk_list:
                        scan_count = write_lcms_chunk_to_mzml(data,
                                                              writer,
                                                              frame_start,
                                                              frame_stop,
                                                              scan_count,
                                                              mode,
                                                              ms2_only,
                                                              exclude_mobility,
                                                              profile_bins,
                                                              mz_encoding,
                                                              intensity_encoding,
                                                              mobility_encoding,
                                                              compression)
                        if gui:
                            sys.stdout.write(get_iso8601_timestamp() +
                                             ':' +
                                             data.source_file.replace('/', '\\') +
                                             ':Progress:100%\n')

    if num_of_spectra != scan_count:
        logging.info(get_iso8601_timestamp() + ':' + 'Updating scan count...')
        update_spectra_count(outdir, outfile, num_of_spectra, scan_count)
    else:
        logging.info(get_iso8601_timestamp() + ':' + 'Renaming mzML file...')
        if os.path.exists(os.path.join(outdir, outfile)):
            os.remove(os.path.join(outdir, outfile))
        os.rename(os.path.splitext(os.path.join(outdir, outfile))[0] + '_tmp.mzML', os.path.join(outdir, outfile))
    logging.info(get_iso8601_timestamp() + ':' + 'Finished writing to .mzML file ' +
                 os.path.join(outdir, outfile) + '...')


def write_maldi_dd_mzml(data, infile, outdir, outfile, mode, ms2_only, exclude_mobility, profile_bins, mz_encoding,
                        intensity_encoding, mobility_encoding, compression, maldi_output_mode, plate_map,
                        barebones_metadata, gui=False):
    """
    Parse and write out spectra to an mzML file from a MALDI-MS(/MS) dried droplet dataset using psims.

    :param data: Object containing raw data information from TDF or TSF file.
    :type data: timsconvert.classes.TimsconvertTdfData | timsconvert.classes.TimsconvertTsfData
    :param infile: Input file path to be used for source file metadata.
    :type infile: str
    :param outdir: Output directory path that was specified from the command line parameters or the original input
        file path if no output directory was specified.
    :type outdir: str
    :param outfile: The original input filename if no output filename was specified.
    :type outfile: str
    :param mode: Mode command line parameter, either "profile", "centroid", or "raw".
    :type mode: str
    :param ms2_only: Whether to include MS1 data in the output files.
    :type ms2_only: bool
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
    :param compression: Compression command line parameter, either "zlib" or "none".
    :type compression: str
    :param maldi_output_mode: Determines whether all spectra from a given .d dataset are written to a single mzML file
        ("combined"), written to individual mzML files ("individual", i.e. one spectrum per file), or grouped to have
        one mzML file per annotation/label/condition ("sample", requires CSV plate_map to be specified).
    :type maldi_output_mode: str
    :param plate_map: Path to the MALDI plate map in CSV format.
    :type plate_map: str
    :param barebones_metadata: If True, omit software and data processing metadata in the resulting mzML files. Used
        for compatibility with downstream analysis software that does not have support for newer CV params or
        UserParams.
    :type barebones_metadata: bool
    :param gui: Whether function is running from GUI or command line. Toggles sys.stdout messages.
    :type gui: bool
    """
    if isinstance(data, TimsconvertBafData):
        frames_key = 'Spectra'
        metadata_key = 'Properties'
    elif isinstance(data, TimsconvertTsfData) or isinstance(data, TimsconvertTdfData):
        frames_key = 'Frames'
        metadata_key = 'GlobalMetadata'
    # All spectra from a given TSF or TDF file are combined into a single mzML file.
    if maldi_output_mode == 'combined':
        # Initialize mzML writer using psims.
        logging.info(get_iso8601_timestamp() + ':' + 'Initializing mzML Writer...')
        writer = MzMLWriter(os.path.splitext(os.path.join(outdir, outfile))[0] + '_tmp.mzML', close=True)

        with writer:
            # Begin mzML with controlled vocabularies (CV).
            logging.info(get_iso8601_timestamp() + ':' + 'Initializing controlled vocabularies...')
            writer.controlled_vocabularies()

            # Start write acquisition, instrument config, processing, etc. to mzML.
            logging.info(get_iso8601_timestamp() + ':' + 'Writing mzML metadata...')
            write_mzml_metadata(data, writer, infile, mode, ms2_only, barebones_metadata)

            logging.info(
                get_iso8601_timestamp() + ':' + 'Writing data to .mzML file ' + os.path.join(outdir, outfile) + '...')
            # Parse chunks of data and write to spectrum element.
            with writer.run(id='run',
                            instrument_configuration='instrument',
                            start_time=data.analysis[metadata_key]['AcquisitionDateTime']):
                scan_count = 0
                # Count number of spectra in run.
                logging.info(get_iso8601_timestamp() + ':' + 'Calculating number of spectra...')
                num_of_spectra = len(data.analysis[frames_key]['Id'].to_list())
                with writer.spectrum_list(count=num_of_spectra):
                    # Parse all MALDI data.
                    num_frames = data.analysis[frames_key].shape[0] + 1
                    # Parse TSF data.
                    if data.analysis[metadata_key]['SchemaType'] == 'TSF':
                        if mode == 'raw':
                            logging.info(get_iso8601_timestamp() + ':' + 'TSF file detected. Only export in profile or '
                                                                         'centroid mode are supported. Defaulting to centroid '
                                                                         'mode.')
                        list_of_scans = parse_maldi_tsf(data,
                                                        1,
                                                        num_frames,
                                                        mode,
                                                        ms2_only,
                                                        profile_bins,
                                                        mz_encoding,
                                                        intensity_encoding)
                    # Parse TDF data.
                    elif data.analysis[metadata_key]['SchemaType'] == 'TDF':
                        list_of_scans = parse_maldi_tdf(data,
                                                        1,
                                                        num_frames,
                                                        mode,
                                                        ms2_only,
                                                        exclude_mobility,
                                                        profile_bins,
                                                        mz_encoding,
                                                        intensity_encoding,
                                                        mobility_encoding)
                    # Write MS1 parent scans.
                    for scan in list_of_scans:
                        if ms2_only and scan.ms_level == 1:
                            pass
                        else:
                            scan_count += 1
                            scan.scan_number = scan_count
                            if scan.ms_level == 1:
                                write_ms1_spectrum(writer,
                                                   data,
                                                   scan,
                                                   mz_encoding,
                                                   intensity_encoding,
                                                   mobility_encoding,
                                                   compression,
                                                   title=os.path.splitext(outfile)[0])
                            elif scan.ms_level == 2:
                                write_ms2_spectrum(writer,
                                                   data,
                                                   scan,
                                                   mz_encoding,
                                                   intensity_encoding,
                                                   mobility_encoding,
                                                   compression,
                                                   title=os.path.splitext(outfile)[0])
                            if gui:
                                sys.stdout.write(get_iso8601_timestamp() +
                                                 ':' +
                                                 data.source_file.replace('/', '\\') +
                                                 ':Progress:' +
                                                 str(round((scan_count / len(list_of_scans)) * 100)) +
                                                 '%\n')

        logging.info(get_iso8601_timestamp() + ':' + 'Updating scan count...')
        update_spectra_count(outdir, outfile, num_of_spectra, scan_count)
        logging.info(
            get_iso8601_timestamp() + ':' + 'Finished writing to .mzML file ' + os.path.join(outdir, outfile) + '...')

    # Each spectrum in a given TSF or TDF file is output as its own individual mzML file.
    elif maldi_output_mode == 'individual' and plate_map != '':
        # Check to make sure plate map is a valid csv file.
        if os.path.exists(plate_map) and os.path.splitext(plate_map)[1] == '.csv':
            # Parse all MALDI data.
            num_frames = data.analysis[frames_key].shape[0] + 1
            # Parse TSF data.
            if data.analysis[metadata_key]['SchemaType'] == 'TSF':
                if mode == 'raw':
                    logging.info(get_iso8601_timestamp() + ':' + 'TSF file detected. Only export in profile or '
                                                                 'centroid mode are supported. Defaulting to centroid '
                                                                 'mode.')
                list_of_scans = parse_maldi_tsf(data,
                                                1,
                                                num_frames,
                                                mode,
                                                ms2_only,
                                                profile_bins,
                                                mz_encoding,
                                                intensity_encoding)
            # Parse TDF data.
            elif data.analysis[metadata_key]['SchemaType'] == 'TDF':
                list_of_scans = parse_maldi_tdf(data,
                                                1,
                                                num_frames,
                                                mode,
                                                ms2_only,
                                                exclude_mobility,
                                                profile_bins,
                                                mz_encoding,
                                                intensity_encoding,
                                                mobility_encoding)

            # Use plate map to determine filename.
            # Names things as sample_position.mzML
            plate_map_dict = parse_maldi_plate_map(plate_map)

            progress_counter = 0
            for scan in list_of_scans:
                progress_counter += 1
                output_filename = os.path.join(outdir,
                                               plate_map_dict[scan.coord] + '_' + scan.coord + '.mzML')

                writer = MzMLWriter(output_filename, close=True)

                with writer:
                    writer.controlled_vocabularies()

                    write_mzml_metadata(data, writer, infile, mode, ms2_only, barebones_metadata)

                    with writer.run(id='run',
                                    instrument_configuration='instrument',
                                    start_time=data.analysis[metadata_key]['AcquisitionDateTime']):
                        scan_count = 1
                        scan.scan_number = scan_count
                        with writer.spectrum_list(count=scan_count):
                            if ms2_only and scan.ms_level == 1:
                                pass
                            else:
                                if scan.ms_level == 1:
                                    write_ms1_spectrum(writer,
                                                       data,
                                                       scan,
                                                       mz_encoding,
                                                       intensity_encoding,
                                                       mobility_encoding,
                                                       compression,
                                                       title=plate_map_dict[scan.coord])
                                elif scan.ms_level == 2:
                                    write_ms2_spectrum(writer,
                                                       data,
                                                       scan,
                                                       mz_encoding,
                                                       intensity_encoding,
                                                       mobility_encoding,
                                                       compression,
                                                       title=plate_map_dict[scan.coord])
                logging.info(get_iso8601_timestamp() + ':' + 'Finished writing to .mzML file ' +
                             os.path.join(outdir, output_filename) + '...')
                if gui:
                    sys.stdout.write(get_iso8601_timestamp() +
                                     ':' +
                                     data.source_file.replace('/', '\\') +
                                     ':Progress:' +
                                     str(round((progress_counter / len(list_of_scans)) * 100)) +
                                     '%\n')

    # Group spectra from a given TSF or TDF file by sample name based on user provided plate map.
    elif maldi_output_mode == 'sample' and plate_map != '':
        # Check to make sure plate map is a valid csv file.
        if os.path.exists(plate_map) and os.path.splitext(plate_map)[1] == '.csv':
            # Parse all MALDI data.
            num_frames = data.analysis[frames_key].shape[0] + 1
            # Parse TSF data.
            if data.analysis[metadata_key]['SchemaType'] == 'TSF':
                if mode == 'raw':
                    logging.info(get_iso8601_timestamp() + ':' + 'TSF file detected. Only export in profile or '
                                                                 'centroid mode are supported. Defaulting to centroid '
                                                                 'mode.')
                list_of_scans = parse_maldi_tsf(data,
                                                1,
                                                num_frames,
                                                mode,
                                                ms2_only,
                                                profile_bins,
                                                mz_encoding,
                                                intensity_encoding)
            # Parse TDF data.
            elif data.analysis[metadata_key]['SchemaType'] == 'TDF':
                list_of_scans = parse_maldi_tdf(data,
                                                1,
                                                num_frames,
                                                mode,
                                                ms2_only,
                                                exclude_mobility,
                                                profile_bins,
                                                mz_encoding,
                                                intensity_encoding,
                                                mobility_encoding)

            # Parse plate map.
            plate_map_dict = parse_maldi_plate_map(plate_map)

            # Get coordinates for each condition replicate.
            conditions = [str(value) for key, value in plate_map_dict.items()]
            conditions = sorted(list(set(conditions)))

            dict_of_scan_lists = {}
            for i in conditions:
                dict_of_scan_lists[i] = []

            for key, value in plate_map_dict.items():
                try:
                    dict_of_scan_lists[value].append(key)
                except KeyError:
                    pass

            progress_counter = 0
            for key, value in dict_of_scan_lists.items():
                if key != 'nan':
                    progress_counter += 1
                    output_filename = os.path.join(outdir, key + '.mzML')

                    writer = MzMLWriter(output_filename, close=True)

                    with writer:
                        writer.controlled_vocabularies()
                        write_mzml_metadata(data, writer, infile, mode, ms2_only, barebones_metadata)
                        with writer.run(id='run',
                                        instrument_configuration='instrument',
                                        start_time=data.analysis[metadata_key]['AcquisitionDateTime']):
                            scan_count = len(value)
                            with writer.spectrum_list(count=scan_count):
                                condition_scans = [i for i in list_of_scans if i.coord in value]
                                scan_count = 0
                                for scan in condition_scans:
                                    if ms2_only and scan.ms_level == 1:
                                        pass
                                    else:
                                        scan_count += 1
                                        scan.scan_number = scan_count
                                        if scan.ms_level == 1:
                                            write_ms1_spectrum(writer,
                                                               data,
                                                               scan,
                                                               mz_encoding,
                                                               intensity_encoding,
                                                               mobility_encoding,
                                                               compression,
                                                               title=key)
                                        elif scan.ms_level == 2:
                                            write_ms2_spectrum(writer,
                                                               data,
                                                               scan,
                                                               mz_encoding,
                                                               intensity_encoding,
                                                               mobility_encoding,
                                                               compression,
                                                               title=key)

                    logging.info(get_iso8601_timestamp() + ':' + 'Finished writing to .mzML file ' +
                                 os.path.join(outdir, outfile) + '...')
                    if gui:
                        sys.stdout.write(get_iso8601_timestamp() +
                                         ':' +
                                         data.source_file.replace('/', '\\') +
                                         ':Progress:' +
                                         str(round((progress_counter / len(dict_of_scan_lists)) * 100)) +
                                         '%\n')


def write_maldi_ims_chunk_to_imzml(data, imzml_file, frame_start, frame_stop, mode, exclude_mobility, profile_bins,
                                   mz_encoding, intensity_encoding, mobility_encoding, diapasef_window=None):
    """
    Parse and write out a group of spectra to an imzML file from a MALDI-MS(/MS) MSI dataset using pyimzML.

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
    """
    # Parse and write TSF data.
    if isinstance(data, TimsconvertTsfData):
        list_of_scans = parse_maldi_tsf(data,
                                        frame_start,
                                        frame_stop,
                                        mode,
                                        False,
                                        profile_bins,
                                        mz_encoding,
                                        intensity_encoding)
        for scans in list_of_scans:
            imzml_file.addSpectrum(scans.mz_array,
                                   scans.intensity_array,
                                   scans.coord)
    # Parse TDF data.
    elif isinstance(data, TimsconvertTdfData):
        if diapasef_window is not None:
            list_of_scans = parse_maldi_tdf_iprm(data,
                                                 frame_start,
                                                 frame_stop,
                                                 mode,
                                                 exclude_mobility,
                                                 profile_bins,
                                                 mz_encoding,
                                                 intensity_encoding,
                                                 mobility_encoding,
                                                 diapasef_window=diapasef_window)
        else:
            list_of_scans = parse_maldi_tdf(data,
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
            for scans in list_of_scans:
                imzml_file.addSpectrum(scans.mz_array,
                                       scans.intensity_array,
                                       scans.coord,
                                       mobilities=scans.mobility_array)
        elif exclude_mobility:
            for scans in list_of_scans:
                imzml_file.addSpectrum(scans.mz_array,
                                       scans.intensity_array,
                                       scans.coord)


def write_maldi_ims_imzml(data, outdir, outfile, mode, exclude_mobility, profile_bins, imzml_mode, mz_encoding,
                          intensity_encoding, mobility_encoding, compression, chunk_size=10, gui=False):
    """
    Parse and write out spectra to an imzML file from a MALDI-MS(/MS) MSI dataset using pyimzML.

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
    :param gui: Whether function is running from GUI or command line. Toggles sys.stdout messages.
    :type gui: bool
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

    logging.info(get_iso8601_timestamp() + ':' + 'Writing to .imzML file ' + os.path.join(outdir, outfile) + '...')
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
                write_maldi_ims_chunk_to_imzml(data,
                                               imzml_file,
                                               frame_start,
                                               frame_stop,
                                               mode,
                                               exclude_mobility,
                                               profile_bins,
                                               mz_encoding,
                                               intensity_encoding,
                                               mobility_encoding)
                if gui:
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
                write_maldi_ims_chunk_to_imzml(data,
                                               imzml_file,
                                               frame_start,
                                               frame_stop,
                                               mode,
                                               exclude_mobility,
                                               profile_bins,
                                               mz_encoding,
                                               intensity_encoding,
                                               mobility_encoding)
                if gui:
                    sys.stdout.write(get_iso8601_timestamp() +
                                     ':' +
                                     data.source_file.replace('/', '\\') +
                                     ':Progress:100%\n')
    logging.info(
        get_iso8601_timestamp() + ':' + 'Finished writing to .imzML file ' + os.path.join(outdir, outfile) + '...')


def write_maldi_ims_iprm_imzml(data, outdir, outfile, mode, exclude_mobility, profile_bins, imzml_mode, mz_encoding,
                               intensity_encoding, mobility_encoding, compression, chunk_size=10, gui=False):
    """
    Parse and write out spectra to an imzML file from a MALDI-MS(/MS) MSI dataset using pyimzML.

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
    :param gui: Whether function is running from GUI or command line. Toggles sys.stdout messages.
    :type gui: bool
    """
    # Set polarity for run in imzML.
    polarity = list(set(data.analysis['Frames']['Polarity'].values.tolist()))
    if len(polarity) == 1 and polarity[0] == '+':
        polarity = 'positive'
    elif len(polarity) == 1 and polarity[0] == '-':
        polarity = 'negative'
    else:
        polarity = None

    # Get compression type object.
    if compression == 'zlib':
        compression_object = ZlibCompression()
    elif compression == 'none':
        compression_object = NoCompression()

    if data.analysis['GlobalMetadata']['SchemaType'] == 'TDF':
        if mode == 'profile':
            exclude_mobility = True
            logging.info(
                get_iso8601_timestamp() + ':' + 'Export of ion mobility data is not supported for profile mode data...')
            logging.info(get_iso8601_timestamp() + ':' + 'Exporting without ion mobility data...')
        msms_mode_id = data.analysis['PropertyDefinitions'][data.analysis['PropertyDefinitions']['PermanentName'] ==
                                                            'Mode_ScanMode'].to_dict(orient='records')[0]['Id']
        properties_dicts = data.analysis['Properties'][data.analysis['Properties']['Property'] ==
                                                       msms_mode_id].to_dict(orient='records')
        if all(i['Value'] == 12 for i in properties_dicts):
            if not data.analysis['DiaFrameMsMsWindows'].empty:
                writers = {}
                for index, row in data.analysis['DiaFrameMsMsWindows'].iterrows():
                    suffix = f"mz{round(row['IsolationMz'])}_ScanNum{round(row['ScanNumBegin'])}-{round(row['ScanNumEnd'])}"
                    if not exclude_mobility:
                        # tuple of ImzMLWriter and diapasef_window
                        writers[suffix] = (ImzMLWriter(os.path.join(outdir,
                                                                    f'{os.path.splitext(outfile)[0]}_{suffix}.imzML'),
                                                       polarity=polarity,
                                                       mode=imzml_mode,
                                                       spec_type=mode,
                                                       mz_dtype=get_encoding_dtype(mz_encoding),
                                                       intensity_dtype=get_encoding_dtype(intensity_encoding),
                                                       mobility_dtype=get_encoding_dtype(mobility_encoding),
                                                       mz_compression=compression_object,
                                                       intensity_compression=compression_object,
                                                       mobility_compression=compression_object,
                                                       include_mobility=True),
                                           row)
                    elif exclude_mobility:
                        writers[suffix] = (ImzMLWriter(os.path.join(outdir,
                                                                    f'{os.path.splitext(outfile)[0]}_{suffix}.imzML'),
                                                       polarity=polarity,
                                                       mode=imzml_mode,
                                                       spec_type=mode,
                                                       mz_dtype=get_encoding_dtype(mz_encoding),
                                                       intensity_dtype=get_encoding_dtype(intensity_encoding),
                                                       mz_compression=compression_object,
                                                       intensity_compression=compression_object,
                                                       include_mobility=False),
                                           row)

                for suffix, value in writers.items():
                    writer = value[0]
                    diapasef_window = value[1]
                    logging.info(get_iso8601_timestamp() + ':' + 'Writing to .imzML file ' + os.path.join(outdir, outfile) + '...')
                    with writer as imzml_file:
                        chunk = 0
                        frames = data.analysis['Frames']['Id'].to_list()
                        while chunk + chunk_size + 1 <= len(frames):
                            chunk_list = []
                            for i, j in zip(frames[chunk:chunk + chunk_size],
                                            frames[chunk + 1: chunk + chunk_size + 1]):
                                chunk_list.append((int(i), int(j)))
                            logging.info(get_iso8601_timestamp() +
                                         ':' +
                                         'Parsing and writing Frame ' +
                                         str(chunk_list[0][0]) +
                                         ' from ' +
                                         data.analysis['GlobalMetadata']['SampleName'] +
                                         '...')
                            for frame_start, frame_stop in chunk_list:
                                write_maldi_ims_chunk_to_imzml(data,
                                                               imzml_file,
                                                               frame_start,
                                                               frame_stop,
                                                               mode,
                                                               exclude_mobility,
                                                               profile_bins,
                                                               mz_encoding,
                                                               intensity_encoding,
                                                               mobility_encoding,
                                                               diapasef_window=diapasef_window)
                                if gui:
                                    sys.stdout.write(get_iso8601_timestamp() +
                                                     ':' +
                                                     data.source_file.replace('/', '\\') +
                                                     ':' +
                                                     suffix +
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
                                write_maldi_ims_chunk_to_imzml(data,
                                                               imzml_file,
                                                               frame_start,
                                                               frame_stop,
                                                               mode,
                                                               exclude_mobility,
                                                               profile_bins,
                                                               mz_encoding,
                                                               intensity_encoding,
                                                               mobility_encoding,
                                                               diapasef_window=diapasef_window)
                                if gui:
                                    sys.stdout.write(get_iso8601_timestamp() +
                                                     ':' +
                                                     data.source_file.replace('/', '\\') +
                                                     ':' +
                                                     suffix +
                                                     ':Progress:100%\n')
                logging.info(
                    get_iso8601_timestamp() + ':' + 'Finished writing to .imzML file ' + os.path.join(outdir, outfile) + '...')


def write_maldi_iprm_ims_chunk_to_mzml(data, mzml_file, frame_start, frame_stop, scan_count, mode, exclude_mobility,
                                       profile_bins, mz_encoding, intensity_encoding, mobility_encoding, compression,
                                       diapasef_window=None):
    # Parse iprm-PASEF TDF data.
    list_of_scans = parse_maldi_tdf_iprm(data,
                                         frame_start,
                                         frame_stop,
                                         mode,
                                         exclude_mobility,
                                         profile_bins,
                                         mz_encoding,
                                         intensity_encoding,
                                         mobility_encoding,
                                         diapasef_window=diapasef_window)
    for scan in list_of_scans:
        if scan.ms_level == 2:
            scan_count += 1
            scan.scan_number = scan_count
            write_ms2_spectrum(mzml_file, data, scan, mz_encoding, intensity_encoding, mobility_encoding, compression)
    return scan_count


def write_maldi_ims_iprm_mzml(data, infile, outdir, outfile, mode, ms2_only, exclude_mobility, profile_bins,
                              mz_encoding, intensity_encoding, mobility_encoding, compression, barebones_metadata,
                              iprm_output_mode, chunk_size=10, gui=False):
    if iprm_output_mode == 'individual':
        if data.analysis['GlobalMetadata']['SchemaType'] == 'TDF':
            if mode == 'profile':
                exclude_mobility = True
                logging.info(get_iso8601_timestamp() + ':' + 'Export of ion mobility data is not supported for profile mode data...')
                logging.info(get_iso8601_timestamp() + ':' + 'Exporting without ion mobility data...')
            msms_mode_id = data.analysis['PropertyDefinitions'][data.analysis['PropertyDefinitions']['PermanentName'] == 'Mode_ScanMode'].to_dict(orient='records')[0]['Id']
            properties_dict = data.analysis['Properties'][data.analysis['Properties']['Property'] == msms_mode_id].to_dict(orient='records')
            if all(i['Value'] == 12 for i in properties_dict):
                if not data.analysis['DiaFrameMsMsWindows'].empty:
                    writers = {}
                    for index, row in data.analysis['DiaFrameMsMsWindows'].iterrows():
                        suffix = f"mz{round(row['IsolationMz'])}_ScanNum{round(row['ScanNumBegin'])}-{round(row['ScanNumEnd'])}"
                        # tuple of MzmlWriter and diapasef_window
                        outfile = f'{os.path.splitext(outfile)[0]}_{suffix}.mzML'
                        writers[suffix] = (MzMLWriter(os.path.join(outdir, outfile), close=True), row)
                    for suffix, value in writers.items():
                        writer = value[0]
                        diapasef_window = value[1]
                        with writer:
                            # Begin mzML with controlled vocabularies (CV).
                            logging.info(get_iso8601_timestamp() + ':' + 'Initializing controlled vocabularies...')
                            writer.controlled_vocabularies()

                            # Start write acquistiion, instrument config, processing, etc. to mzML.
                            logging.info(get_iso8601_timestamp() + ':' + 'Writing mzML metadata...')
                            write_mzml_metadata(data, writer, infile, mode, ms2_only, barebones_metadata)

                            logging.info(get_iso8601_timestamp() + ':' + 'Writing to .mzML file ' + os.path.join(outdir, outfile) + '...')
                            # Parse chunks of data and write to spectrum element.
                            with writer.run(id='run',
                                            instrument_configuration='instrument',
                                            start_time=data.analysis['GlobalMetadata']['AcquisitionDateTime']):
                                scan_count = 0
                                # Count number of spectra in run by counting number of rows in Frames table.
                                logging.info(get_iso8601_timestamp() + ':' + 'Calculating number of spectra...')
                                num_of_spectra = len(data.analysis['Frames']['Id'].to_list())
                                with writer.spectrum_list(count=num_of_spectra):
                                    chunk = 0
                                    frames = data.analysis['Frames']['Id'].to_list()
                                    # Write data in chunks of chunk_size.
                                    while chunk + chunk_size + 1 <= len(frames):
                                        chunk_list = []
                                        for i, j in zip(frames[chunk:chunk + chunk_size],
                                                        frames[chunk + 1: chunk + chunk_size + 1]):
                                            chunk_list.append((int(i), int(j)))
                                        logging.info(get_iso8601_timestamp() +
                                                     ':' +
                                                     'Parsing and writing Frame ' +
                                                     str(chunk_list[0][0]) +
                                                     ' from ' +
                                                     data.analysis['GlobalMetadata']['SampleName'] +
                                                     '...')
                                        for frame_start, frame_stop in chunk_list:
                                            scan_count = write_maldi_iprm_ims_chunk_to_mzml(data,
                                                                                            writer,
                                                                                            frame_start,
                                                                                            frame_stop,
                                                                                            scan_count,
                                                                                            mode,
                                                                                            exclude_mobility,
                                                                                            profile_bins,
                                                                                            mz_encoding,
                                                                                            intensity_encoding,
                                                                                            mobility_encoding,
                                                                                            compression,
                                                                                            diapasef_window)
                                            if gui:
                                                sys.stdout.write(get_iso8601_timestamp() +
                                                                 ':' +
                                                                 data.source_file.replace('/', '\\') +
                                                                 ':' +
                                                                 suffix +
                                                                 ':Progress:' +
                                                                 str(round((frame_start / data.analysis['Frames'].shape[0]) * 100)) +
                                                                 '%\n')
                                        chunk += chunk_size
                                    # Last chunk may be smaller than chunk_size.
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
                                            scan_count = write_maldi_iprm_ims_chunk_to_mzml(data,
                                                                                            writer,
                                                                                            frame_start,
                                                                                            frame_stop,
                                                                                            scan_count,
                                                                                            mode,
                                                                                            exclude_mobility,
                                                                                            profile_bins,
                                                                                            mz_encoding,
                                                                                            intensity_encoding,
                                                                                            mobility_encoding,
                                                                                            compression,
                                                                                            diapasef_window)
                                            if gui:
                                                sys.stdout.write(get_iso8601_timestamp() +
                                                                 ':' +
                                                                 data.source_file.replace('/', '\\') +
                                                                 ':' +
                                                                 suffix +
                                                                 ':Progress:100%\n')
                        logging.info(get_iso8601_timestamp() + ':' + 'Finished writing to .mzML file ' +
                                     os.path.join(outdir, outfile) + '...')
    elif iprm_output_mode == 'combined':
        if data.analysis['GlobalMetadata']['SchemaType'] == 'TDF':
            if mode == 'profile':
                exclude_mobility = True
                logging.info(
                    get_iso8601_timestamp() + ':' + 'Export of ion mobility data is not supported for profile mode data...')
                logging.info(get_iso8601_timestamp() + ':' + 'Exporting without ion mobility data...')
            msms_mode_id = data.analysis['PropertyDefinitions'][
                data.analysis['PropertyDefinitions']['PermanentName'] == 'Mode_ScanMode'].to_dict(orient='records')[0][
                'Id']
            properties_dict = data.analysis['Properties'][
                data.analysis['Properties']['Property'] == msms_mode_id].to_dict(orient='records')
            if all(i['Value'] == 12 for i in properties_dict):
                if not data.analysis['DiaFrameMsMsWindows'].empty:
                    # Initialize psims mzML writer.
                    writer = MzMLWriter(os.path.join(outdir, f'{os.path.splitext(outfile)[0]}.mzML'), close=True)
                    with writer:
                        # Begin mzML with controlled vocabularies (CV).
                        logging.info(get_iso8601_timestamp() + ':' + 'Initializing controlled vocabularies...')
                        writer.controlled_vocabularies()

                        # Start write acquistiion, instrument config, processing, etc. to mzML.
                        logging.info(get_iso8601_timestamp() + ':' + 'Writing mzML metadata...')
                        write_mzml_metadata(data, writer, infile, mode, ms2_only, barebones_metadata)

                        logging.info(
                            get_iso8601_timestamp() + ':' + 'Writing to .mzML file ' + os.path.join(outdir, outfile) + '...')
                        # Parse chunks of data and write to spectrum element.
                        with writer.run(id='run',
                                        instrument_configuration='instrument',
                                        start_time=data.analysis['GlobalMetadata']['AcquisitionDateTime']):
                            scan_count = 0
                            # Count number of spectra in run by counting number of rows in Frames table.
                            logging.info(get_iso8601_timestamp() + ':' + 'Calculating number of spectra...')
                            num_of_spectra = len(data.analysis['Frames']['Id'].to_list())
                            with writer.spectrum_list(count=num_of_spectra):
                                chunk = 0
                                frames = data.analysis['Frames']['Id'].to_list()
                                # Write data in chunks of chunk_size.
                                while chunk + chunk_size + 1 <= len(frames):
                                    chunk_list = []
                                    for i, j in zip(frames[chunk:chunk + chunk_size],
                                                    frames[chunk + 1: chunk + chunk_size + 1]):
                                        chunk_list.append((int(i), int(j)))
                                    logging.info(get_iso8601_timestamp() +
                                                 ':' +
                                                 'Parsing and writing Frame ' +
                                                 str(chunk_list[0][0]) +
                                                 ' from ' +
                                                 data.analysis['GlobalMetadata']['SampleName'] +
                                                 '...')
                                    for frame_start, frame_stop in chunk_list:
                                        for index, diapasef_window in data.analysis['DiaFrameMsMsWindows'].iterrows():
                                            scan_count = write_maldi_iprm_ims_chunk_to_mzml(data,
                                                                                            writer,
                                                                                            frame_start,
                                                                                            frame_stop,
                                                                                            scan_count,
                                                                                            mode,
                                                                                            exclude_mobility,
                                                                                            profile_bins,
                                                                                            mz_encoding,
                                                                                            intensity_encoding,
                                                                                            mobility_encoding,
                                                                                            compression,
                                                                                            diapasef_window)
                                            if gui:
                                                sys.stdout.write(get_iso8601_timestamp() +
                                                                 ':' +
                                                                 data.source_file.replace('/', '\\') +
                                                                 ':Progress:' +
                                                                 str(round((frame_start / data.analysis['Frames'].shape[0]) * 100)) +
                                                                 '%\n')
                                    chunk += chunk_size
                                # Last chunk may be smaller than chunk_size
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
                                        for index, diapasef_window in data.analysis['DiaFrameMsMsWindows'].iterrows():
                                            scan_count = write_maldi_iprm_ims_chunk_to_mzml(data,
                                                                                            writer,
                                                                                            frame_start,
                                                                                            frame_stop,
                                                                                            scan_count,
                                                                                            mode,
                                                                                            exclude_mobility,
                                                                                            profile_bins,
                                                                                            mz_encoding,
                                                                                            intensity_encoding,
                                                                                            mobility_encoding,
                                                                                            compression,
                                                                                            diapasef_window)
                                            if gui:
                                                sys.stdout.write(get_iso8601_timestamp() +
                                                                 ':' +
                                                                 data.source_file.replace('/', '\\') +
                                                                 ':Progress:100%\n')
                    logging.info(get_iso8601_timestamp() + ':' + 'Finished writing to .mzML file ' +
                                 os.path.join(outdir, outfile) + '...')


def write_maldi_ims_iprm_mgf(data, outdir, outfile, mode, exclude_mobility, profile_bins, mz_encoding,
                             intensity_encoding, mobility_encoding, iprm_output_mode, chunk_size=10, gui=False):
    # TODO: set exclude mobility to True for mgf output
    if iprm_output_mode == 'individual':
        if data.analysis['GlobalMetadata']['SchemaType'] == 'TDF':
            if mode == 'profile':
                exclude_mobility = True
                logging.info(
                    get_iso8601_timestamp() + ':' + 'Export of ion mobility data is not supported for profile mode data...')
                logging.info(get_iso8601_timestamp() + ':' + 'Exporting without ion mobility data...')
            msms_mode_id = data.analysis['PropertyDefinitions'][
                data.analysis['PropertyDefinitions']['PermanentName'] == 'Mode_ScanMode'].to_dict(orient='records')[0][
                'Id']
            properties_dict = data.analysis['Properties'][
                data.analysis['Properties']['Property'] == msms_mode_id].to_dict(orient='records')
            if all(i['Value'] == 12 for i in properties_dict):
                if not data.analysis['DiaFrameMsMsWindows'].empty:
                    writers = {}
                    for index, row in data.analysis['DiaFrameMsMsWindows'].iterrows():
                        suffix = f"mz{round(row['IsolationMz'])}_ScanNum{round(row['ScanNumBegin'])}-{round(row['ScanNumEnd'])}"
                        # tuple of mgf filename and diapasef_window
                        outfile = f'{os.path.splitext(outfile)[0]}_{suffix}.mgf'
                        writers[suffix] = (os.path.join(outdir, outfile), row)
                    for suffix, value in writers.items():
                        mgf_filename = value[0]
                        diapasef_window = value[1]
                        scan_count = 0
                        chunk = 0
                        frames = data.analysis['Frames']['Id'].to_list()
                        # Write data in chunks of chunk_size.
                        while chunk + chunk_size + 1 <= len(frames):
                            chunk_list = []
                            for i, j in zip(frames[chunk:chunk + chunk_size],
                                            frames[chunk + 1: chunk + chunk_size + 1]):
                                chunk_list.append((int(i), int(j)))
                            logging.info(get_iso8601_timestamp() +
                                         ':' +
                                         'Parsing and writing Frame ' +
                                         str(chunk_list[0][0]) +
                                         ' from ' +
                                         data.analysis['GlobalMetadata']['SampleName'] +
                                         '...')
                            for frame_start, frame_stop in chunk_list:
                                # Parse iprm-PASEF TDF data.
                                list_of_scans = parse_maldi_tdf_iprm(data,
                                                                     frame_start,
                                                                     frame_stop,
                                                                     mode,
                                                                     exclude_mobility,
                                                                     profile_bins,
                                                                     mz_encoding,
                                                                     intensity_encoding,
                                                                     mobility_encoding,
                                                                     diapasef_window=diapasef_window)
                                if scan_count == 0:
                                    mgf_file_mode = 'w'
                                else:
                                    mgf_file_mode = 'a'
                                scan_count += 1
                                ms2_dict_list = [{'m/z array': copy.deepcopy(scan.mz_array),
                                                  'intensity array': copy.deepcopy(scan.intensity_array),
                                                  'params': {'FEATURE_ID': scan_count,
                                                             'PEPMASS': scan.selected_ion_mz,
                                                             'ION_MOBILITY': scan.selected_ion_mobility,
                                                             'MSLEVEL': 2}}
                                                 for scan in list_of_scans]
                                mgf.write(ms2_dict_list, output=mgf_filename, file_mode=mgf_file_mode)
                                if gui:
                                    sys.stdout.write(get_iso8601_timestamp() +
                                                     ':' +
                                                     data.source_file.replace('/', '\\') +
                                                     ':' +
                                                     suffix +
                                                     ':Progress:' +
                                                     str(round((frame_start / data.analysis['Frames'].shape[0]) * 100)) +
                                                     '%\n')
                            chunk += chunk_size
                        # Last chunk may be smaller than chunk_size.
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
                                # Parse iprm-PASEF TDF data.
                                list_of_scans = parse_maldi_tdf_iprm(data,
                                                                     frame_start,
                                                                     frame_stop,
                                                                     mode,
                                                                     exclude_mobility,
                                                                     profile_bins,
                                                                     mz_encoding,
                                                                     intensity_encoding,
                                                                     mobility_encoding,
                                                                     diapasef_window=diapasef_window)
                                if scan_count == 0:
                                    mgf_file_mode = 'w'
                                else:
                                    mgf_file_mode = 'a'
                                scan_count += 1
                                ms2_dict_list = [{'m/z array': copy.deepcopy(scan.mz_array),
                                                  'intensity array': copy.deepcopy(scan.intensity_array),
                                                  'params': {'FEATURE_ID': scan_count,
                                                             'PEPMASS': scan.selected_ion_mz,
                                                             'ION_MOBILITY': scan.selected_ion_mobility,
                                                             'MSLEVEL': 2}}
                                                 for scan in list_of_scans]
                                mgf.write(ms2_dict_list, output=mgf_filename, file_mode=mgf_file_mode)
                                if gui:
                                    sys.stdout.write(get_iso8601_timestamp() +
                                                     ':' +
                                                     data.source_file.replace('/', '\\') +
                                                     ':' +
                                                     suffix +
                                                     ':Progress:100%\n')
                        logging.info(get_iso8601_timestamp() + ':' + 'Finished writing to .mgf file ' +
                                     os.path.join(outdir, outfile) + '...')
    elif iprm_output_mode == 'combined':
        if data.analysis['GlobalMetadata']['SchemaType'] == 'TDF':
            if mode == 'profile':
                exclude_mobility = True
                logging.info(
                    get_iso8601_timestamp() + ':' + 'Export of ion mobility data is not supported for profile mode data...')
                logging.info(get_iso8601_timestamp() + ':' + 'Exporting without ion mobility data...')
            msms_mode_id = data.analysis['PropertyDefinitions'][
                data.analysis['PropertyDefinitions']['PermanentName'] == 'Mode_ScanMode'].to_dict(orient='records')[0][
                'Id']
            properties_dict = data.analysis['Properties'][
                data.analysis['Properties']['Property'] == msms_mode_id].to_dict(orient='records')
            if all(i['Value'] == 12 for i in properties_dict):
                if not data.analysis['DiaFrameMsMsWindows'].empty:
                    mgf_filename = os.path.join(outdir, f'{os.path.splitext(outfile)[0]}.mgf')
                    scan_count = 0
                    chunk = 0
                    frames = data.analysis['Frames']['Id'].to_list()
                    # Write data in chunks of chunk_size.
                    while chunk + chunk_size + 1 <= len(frames):
                        chunk_list = []
                        for i, j in zip(frames[chunk:chunk + chunk_size],
                                        frames[chunk + 1: chunk + chunk_size + 1]):
                            chunk_list.append((int(i), int(j)))
                        logging.info(get_iso8601_timestamp() +
                                     ':' +
                                     'Parsing and writing Frame ' +
                                     str(chunk_list[0][0]) +
                                     ' from ' +
                                     data.analysis['GlobalMetadata']['SampleName'] +
                                     '...')
                        for frame_start, frame_stop in chunk_list:
                            for index, diapasef_window in data.analysis['DiaFrameMsMsWindows'].iterrows():
                                # Parse iprm-PASEF TDF data.
                                list_of_scans = parse_maldi_tdf_iprm(data,
                                                                     frame_start,
                                                                     frame_stop,
                                                                     mode,
                                                                     exclude_mobility,
                                                                     profile_bins,
                                                                     mz_encoding,
                                                                     intensity_encoding,
                                                                     mobility_encoding,
                                                                     diapasef_window=diapasef_window)
                                if scan_count == 0:
                                    mgf_file_mode = 'w'
                                else:
                                    mgf_file_mode = 'a'
                                scan_count += 1
                                ms2_dict_list = [{'m/z array': copy.deepcopy(scan.mz_array),
                                                  'intensity array': copy.deepcopy(scan.intensity_array),
                                                  'params': {'FEATURE_ID': scan_count,
                                                             'PEPMASS': scan.selected_ion_mz,
                                                             'ION_MOBILITY': scan.selected_ion_mobility,
                                                             'MSLEVEL': 2}}
                                                 for scan in list_of_scans]
                                mgf.write(ms2_dict_list, output=mgf_filename, file_mode=mgf_file_mode)
                                if gui:
                                    sys.stdout.write(get_iso8601_timestamp() +
                                                     ':' +
                                                     data.source_file.replace('/', '\\') +
                                                     ':Progress:' +
                                                     str(round((frame_start / data.analysis['Frames'].shape[0]) * 100)) +
                                                     '%\n')
                        chunk += chunk_size
                    # Last chunk may be smaller than chunk_size.
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
                            for index, diapasef_window in data.analysis['DiaFrameMsMsWindows'].iterrows():
                                # Parse iprm-PASEF TDF data.
                                list_of_scans = parse_maldi_tdf_iprm(data,
                                                                     frame_start,
                                                                     frame_stop,
                                                                     mode,
                                                                     exclude_mobility,
                                                                     profile_bins,
                                                                     mz_encoding,
                                                                     intensity_encoding,
                                                                     mobility_encoding,
                                                                     diapasef_window=diapasef_window)
                                if scan_count == 0:
                                    mgf_file_mode = 'w'
                                else:
                                    mgf_file_mode = 'a'
                                scan_count += 1
                                ms2_dict_list = [{'m/z array': copy.deepcopy(scan.mz_array),
                                                  'intensity array': copy.deepcopy(scan.intensity_array),
                                                  'params': {'FEATURE_ID': scan_count,
                                                             'PEPMASS': scan.selected_ion_mz,
                                                             'ION_MOBILITY': scan.selected_ion_mobility,
                                                             'MSLEVEL': 2}}
                                                 for scan in list_of_scans]
                                mgf.write(ms2_dict_list, output=mgf_filename, file_mode=mgf_file_mode)
                                if gui:
                                    sys.stdout.write(get_iso8601_timestamp() +
                                                     ':' +
                                                     data.source_file.replace('/', '\\') +
                                                     ':Progress:100%\n')
                    logging.info(get_iso8601_timestamp() + ':' + 'Finished writing to .mgf file ' +
                                 os.path.join(outdir, outfile) + '...')
