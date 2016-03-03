##########################################################################
# NSAp - Copyright (C) CEA, 2015 - 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

# System import
import os
import gzip
import shutil

# Clindmri import
from clindmri.extensions.connectomist.exceptions import (
    ConnectomistBadFileError)
from clindmri.extensions.connectomist import PtkWrapper


def ptk_nifti_to_gis(nifti, gis, nb_tries=10):
    """
    Function that wraps the PtkNifti2GisConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    nifti: str
        path to the input Nifti file to be converted.
    gis: str
        path without extension to the 3 output GIS files.
    nb_tries: int
        number of times to try the conversion. It often fails
        when using parallel processing.

    Returns
    -------
    gis: str
        path to the generated Gis file.

    Raises
    ------
    ConnectomistError: If call to PtkNifti2GisConverter failed nb_tries times.
    """
    # Check input existence
    if not os.path.isfile(nifti):
        raise ConnectomistBadFileError(nifti)

    # Add extension if there is none
    if not gis.endswith("ima"):
        gis += ".ima"

    # Call command line tool
    cmd = ["PtkNifti2GisConverter", "-i", nifti, "-o", gis,
           "-verbose", "False", "-verbosePluginLoading", "False"]
    files_to_check = [gis]
    ptkprocess = PtkWrapper(cmd)
    ptkprocess(files_to_check, nb_tries=nb_tries)

    return gis


def gz_compress(file_to_compress, clean=True):
    """
    Compress a file with gzip, the output path is the same but with
    ".gz" extension.

    The function raises ConnectomistBadFileError if the input file does not
    exist or if the output compressed file is not created.

    Parameters
    ----------
    file_to_compress: str
        the file to compress.
    clean: bool (optional, default True)
        If 'clean' is True, the input file is deleted, to keep only the
        compressed version.

    Returns
    -------
    gz_file: str
        the compressed output file.
    """
    # Check if the input file exists
    if not os.path.isfile(file_to_compress):
        raise ConnectomistBadFileError(file_to_compress)

    # Zip the input file
    gz_file = file_to_compress + ".gz"
    with open(file_to_compress, 'rb') as f_in, gzip.open(gz_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    if clean:
        os.remove(file_to_compress)
    if not os.path.isfile(gz_file):
        raise ConnectomistBadFileError(gz_file)

    return gz_file


def ptk_gis_to_nifti(gis, nifti, nb_tries=10):
    """
    Function that wraps the PtkGis2NiftiConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    gis: str
        path to the Gis .ima file.
    nifti: str
        path to the output Nifti file.
    nb_tries: int (optional, default 10)
        number of times to try the conversion. It often fails
        when using parallel processing.

    Returns
    -------
    nifti: str
        path to output file.

    Raises
    ------
    ConnectomistError: If call to PtkGis2NiftiConverter failed nb_tries times.
    """
    # Check input existence
    if not os.path.isfile(gis):
        raise ConnectomistBadFileError(gis)

    # The command line tool does not handle .gz properly
    compress_to_gz = False
    if nifti.endswith(".gz"):
        nifti = nifti[:-3]  # remove ".gz"
        compress_to_gz = True

    # Add .nii extension if not the case
    if not nifti.endswith(".nii"):
        nifti += ".nii"

    # Call command line tool:
    # it creates a Nifti + a .minf file (metainformation)
    cmd = ["PtkGis2NiftiConverter", "-i", gis, "-o", nifti,
           "-verbose", "False", "-verbosePluginLoading", "False"]
    files_to_check = [nifti]
    ptkprocess = PtkWrapper(cmd)
    ptkprocess(files_to_check, nb_tries=nb_tries)

    # Compress to nifti if requested
    if compress_to_gz:
        nifti = gz_compress(nifti)

    return nifti


def ptk_concatenate_volumes(path_inputs, path_output, axis="t", nb_tries=10):
    """
    Function that wraps the PtkCat command line tool from Connectomist, with
    only the basic arguments. It allows concatenating volumes. In particular to
    concatenante T2 and DW volumes in one file at the end of the preprocessing.

    Parameters
    ----------
    path_inputs: list of str
        paths to input volumes.
    path_output: str
        path to output volume.
    axis: str (optional, default 't')
        axis along which the concatenation is done.
    nb_tries: int (optional, default 10)
        number of times to try the concatenation.
        The reason why it sometimes fails is not known.

    Returns
    -------
    path_output: str
        path to output gis volume.

    Raises
    ------
    ConnectomistError: If call to PtkCat failed.
    """
    # Check input existence
    for path in path_inputs:
        if not os.path.isfile(path):
            raise ConnectomistBadFileError(path)

    # Add extension if there is none
    if not path_output.endswith("ima"):
        path_output += ".ima"

    # Run command line tool
    cmd = ["PtkCat", "-i"] + path_inputs + [
        "-o", path_output, "-t", axis, "-verbose", "False",
        "-verbosePluginLoading", "False"]
    files_to_check = [path_output]
    ptkprocess = PtkWrapper(cmd)
    ptkprocess(files_to_check, nb_tries=nb_tries)

    return path_output


def ptk_split_t2_and_diffusion(t2_dw_input, t2_output, dw_output, nb_tries=10):
    """
    Function meant to split a Gis file containing a T2 volume (first volume)
    and diffusion-weigthed volumes (the other volumes) in 2 Gis files.
    The separation is done using the PtkSubVolume command line tool from
    Connectomist.

    Parameters
    ----------
    t2_dw_input: str
        path to input volume.
    t2_output: str
        path to output T2 volume.
    dw_output: str
        path to output diffusion-weighted volumes.
    nb_tries: int
        number of times to try the separation. The reason why it sometimes
        fails is not known.

    Returns
    -------
    t2_output: str
        the t2 diffusion Gis volume.
    dw_output: str
        the diffusion weighted Gis volume.

    Raises
    ------
    ConnectomistError: If call to PtkSubVolume failed.
    """
    # Check input existence
    if not os.path.isfile(t2_dw_input):
        raise ConnectomistBadFileError(t2_dw_input)

    # Check that input is a Gis file
    if not t2_dw_input.endswith(".ima"):
        raise ConnectomistBadFileError(t2_dw_input)

    # Set the Gis extension to output file
    if not t2_output.endswith(".ima"):
        t2_output += ".ima"
    if not dw_output.endswith(".ima"):
        dw_output += ".ima"

    # Step 1 - extract the T2 (nodif volume), assuming only one volume with
    # bvalue=0
    cmd_t2 = ["PtkSubVolume", "-i", t2_dw_input, "-o", t2_output, "-tIndices",
              "0", "-verbose", "false", "-verbosePluginLoading", "false"]
    files_to_check = [t2_output]
    ptkprocess = PtkWrapper(cmd_t2)
    ptkprocess(files_to_check, nb_tries=nb_tries)

    # Step 2 - extract the diffusion volumes
    # assuming only all volumes except the first one are diffusion-weighted
    cmd_dw = ["PtkSubVolume", "-i", t2_dw_input, "-o", dw_output, "-t", "1",
              "-verbose", "false", "-verbosePluginLoading", "false"]
    files_to_check = [dw_output]
    ptkprocess = PtkWrapper(cmd_dw)
    ptkprocess(files_to_check, nb_tries=nb_tries)

    return t2_output, dw_output
