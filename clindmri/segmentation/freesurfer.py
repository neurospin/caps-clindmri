#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013-2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import glob
import nibabel
import numpy
import json
import csv
import shutil
from nibabel import freesurfer

# Clindmri import
from clindmri.extensions.fsl.exceptions import FSLResultError
from clindmri.extensions.freesurfer.exceptions import FreeSurferRuntimeError
from clindmri.extensions.freesurfer.wrappers import FSWrapper
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.extensions.freesurfer.reader import TriSurface
from clindmri.extensions.freesurfer.reader import apply_affine_on_mesh
from clindmri.extensions.freesurfer.reader import tkregister_translation
from clindmri.registration.fsl import flirt
from clindmri.extensions.fsl import flirt2aff
from clindmri.registration.utils import extract_image


"""
create a registration matrix between the conformed space (orig.mgz) and the
native anatomical (rawavg.mgz)

tkregister2 --mov rawavg.mgz --targ orig.mgz --reg register.native.dat
--noedit --regheader

map the surface to the native space:

mri_surf2surf --sval-xyz pial --reg register.native.dat rawavg.mgz --tval
lh.pial.native --tval-xyz --hemi lh --s subjectname
"""


def recon_all(fsdir, anatfile, sid, output_directory=None,
              fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Performs all the FreeSurfer cortical reconstruction process.

    Processing stages:

    * Motion Correction and Conform
    * NU (Non-Uniform intensity normalization)
    * Talairach transform computation
    * Intensity Normalization 1
    * Skull Strip
    * EM Register (linear volumetric registration)
    * CA Intensity Normalization
    * CA Non-linear Volumetric Registration
    * Remove Neck
    * LTA with Skull
    * CA Label (Volumetric Labeling, ie Aseg) and Statistics
    * Intensity Normalization 2 (start here for control points)
    * White matter segmentation
    * Edit WM With ASeg
    * Fill (start here for wm edits)
    * Tessellation (begins per-hemisphere operations)
    * Smooth1
    * Inflate1
    * QSphere
    * Automatic Topology Fixer
    * Final Surfs (start here for brain edits for pial surf)
    * Smooth2
    * Inflate2
    * Spherical Mapping
    * Spherical Registration
    * Spherical Registration, Contralateral hemisphere
    * Map average curvature to subject
    * Cortical Parcellation - Desikan_Killiany and Christophe (Labeling)
    * Cortical Parcellation Statistics
    * Cortical Ribbon Mask
    * Cortical Parcellation mapping to Aseg

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="anatfile" type="File" desc="The input anatomical image
            to be segmented with freesurfer."/>
        <input name="sid" type="Str" description="The current subject
            identifier."/>
        <input name="output_directory" type="Directory" description="The
            freesurfer runtime folder."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="subjfsdir" type="Directory" description="Path to the
            resulting freesurfer segmentation."/>
    </unit>
    """
    # Create the fs output directory if necessary
    if not os.path.isdir(fsdir):
        os.makedirs(fsdir)

    # Call freesurfer
    cmd = ["recon-all", "-all", "-subjid", sid, "-i", anatfile, "-sd", fsdir]
    recon = FSWrapper(cmd, shfile=fsconfig)
    recon()
    if recon.exitcode != 0:
        raise FreeSurferRuntimeError(
            recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr + recon.stdout)
    subjfsdir = os.path.join(fsdir, sid)

    return subjfsdir


def aparcstats2table(fsdir, output_directory=None,
                     fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Generate text/ascii tables of freesurfer parcellation stats data
    '?h.aparc.stats'. This can then be easily imported into a spreadsheet
    and/or stats program.

    The labels are located here: $FREESURFER_HOME/FreeSurferColorLUT.txt

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="output_directory" type="Directory" description="The
            statistical destination folder."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="statfiles" type="List_File" description="The freesurfer
            summary stats."/>
    </unit>
    """
    # Parameter that will contain the output stats
    statfiles = []

    # Fist find all the subjects with a stat dir
    statdirs = glob.glob(os.path.join(fsdir, "*", "stats"))
    subjects = [item.lstrip(os.sep).split("/")[-2] for item in statdirs]
    if output_directory is not None:
        path = os.path.join(output_directory, "subjects.json")
        with open(path, "w") as open_file:
            json.dump(subjects, open_file, indent=4)

    # Save the FreeSurfer current working directory and set the new one
    fscwd = None
    if "SUBJECTS_DIR" in os.environ:
        fscwd = os.environ["SUBJECTS_DIR"]
    os.environ["SUBJECTS_DIR"] = fsdir

    # Create the output stat directory
    fsoutdir = os.path.join(fsdir, "stats")
    if not os.path.isdir(fsoutdir):
        os.makedirs(fsoutdir)

    # Call freesurfer
    for hemi in ["lh", "rh"]:
        for meas in ["area", "volume", "thickness", "thicknessstd",
                     "meancurv", "gauscurv", "foldind", "curvind"]:

            statfile = os.path.join(
                fsoutdir, "aparc_stats_{0}_{1}.csv".format(hemi, meas))
            statfiles.append(statfile)
            cmd = ["aparcstats2table", "--subjects"] + subjects + [
                "--hemi", hemi, "--meas", meas, "--tablefile", statfile,
                "--delimiter", "comma", "--parcid-only"]

            recon = FSWrapper(cmd, shfile=fsconfig)
            recon()
            if recon.exitcode != 0:
                raise FreeSurferRuntimeError(
                    recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr +
                    recon.stdout)

    # Restore the FreeSurfer working directory
    if fscwd is not None:
        os.environ["SUBJECTS_DIR"] = fscwd

    return statfiles


def asegstats2table(fsdir, output_directory=None,
                    fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Generate text/ascii tables of freesurfer parcellation stats data
    'aseg.stats'. This can then be easily imported into a spreadsheet
    and/or stats program.

    The labels are located here: $FREESURFER_HOME/FreeSurferColorLUT.txt

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="output_directory" type="Directory" description="The
            statistical destination folder."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="statfiles" type="List_File" description="The freesurfer
            summary stats."/>
    </unit>
    """
    # Parameter that will contain the output stats
    statfiles = []

    # Fist find all the subjects with a stat dir
    statdirs = glob.glob(os.path.join(fsdir, "*", "stats"))
    subjects = [item.lstrip(os.sep).split("/")[-2] for item in statdirs]
    if output_directory is not None:
        path = os.path.join(output_directory, "subjects.json")
        with open(path, "w") as open_file:
            json.dump(subjects, open_file, indent=4)

    # Save the FreeSurfer current working directory and set the new one
    fscwd = None
    if "SUBJECTS_DIR" in os.environ:
        fscwd = os.environ["SUBJECTS_DIR"]
    os.environ["SUBJECTS_DIR"] = fsdir

    # Create the output stat directory
    fsoutdir = os.path.join(fsdir, "stats")
    if not os.path.isdir(fsoutdir):
        os.makedirs(fsoutdir)

    # Call freesurfer
    statfile = os.path.join(fsoutdir, "aseg_stats_volume.csv")
    statfiles.append(statfile)
    cmd = ["asegstats2table", "--subjects"] + subjects + [
        "--meas", "volume", "--tablefile", statfile, "--delimiter", "comma"]
    recon = FSWrapper(cmd, shfile=fsconfig)
    recon()
    if recon.exitcode != 0:
        raise FreeSurferRuntimeError(
            recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr +
            recon.stdout)

    # Restore the FreeSurfer working directory
    if fscwd is not None:
        os.environ["SUBJECTS_DIR"] = fscwd

    return statfiles


def textures2table(regex, ico_order, fsdir, output_directory=None,
                   keep_individual_textures=False, save_mode="numpy",
                   fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Generate text/ascii tables of freesurfer textures data. This can then
    be easily imported into a spreadsheet and/or stats program. Note that
    all the subject texture vertices need to be resampled in a common space.

    Parameters
    ----------
    regex: str (mandatory)
        a regular expression used to locate the files to be converted from
        the 'fsdir' directory.
    ico_order: int (mandatory)
        icosahedron order in [0, 7] that will be used to generate the cortical
        surface texture at a specific tessalation (the corresponding cortical
        surface can be resampled using the
        'clindmri.segmentation.freesurfer.resample_cortical_surface' function).
    fsdir: str (mandatory)
        FreeSurfer subjects directory 'SUBJECTS_DIR'.
    output_directory: str (optional, default None)
        a folder where some information about the processing could
        be written.
    keep_individual_textures: bool (optional, default False)
        if True, keep the individual resampled subject textures on disk.
    save_mode: str (optional, default 'numpy')
        save result in 'csv' or in 'numpy' or 'all'. In CSV format we keep
        only 4 digits and in Numpy format we save several single arrays into
        a single file.
    fsconfig: str (optional)
        the FreeSurfer '.sh' config file.

    Returns
    -------
    textures_files: list of str
        a list of file containing the selected subjects summary texture values.
    """
    # Check input parameters
    if save_mode not in ["numpy", "csv", "all"]:
        raise ValueError("'{0}' is not a valid save option must be "
                         "in ['numpy', 'csv', 'all']".format(save_mode))

    # Get the requested subject textures from the regex
    textures = glob.glob(os.path.join(fsdir, regex))
    if output_directory is not None:
        path = os.path.join(output_directory, "textrues.json")
        with open(path, "w") as open_file:
            json.dump(textures, open_file, indent=4)

    # Resample each texture file
    basename = os.path.basename(regex)
    fsoutdir = os.path.join(fsdir, "textures")
    surfacesdir = os.path.join(fsoutdir, basename)
    if not os.path.isdir(surfacesdir):
        os.makedirs(surfacesdir)
    textures_map = {}
    for texturefile in textures:

        # Get the subject id
        sid = texturefile.replace(fsdir, "")
        sid = sid.lstrip(os.sep).split(os.sep)[0]

        # Get the hemisphere
        hemi = basename.split(".")[0]

        # Create the destination resamples texture file
        resampled_texturefile = os.path.join(surfacesdir, "{0}_{1}.mgz".format(
            sid, basename))

        # Reasmple the surface
        mri_surf2surf(hemi, texturefile, resampled_texturefile,
                      ico_order=ico_order, fsdir=fsdir, sid=sid,
                      fsconfig=fsconfig)

        # Check that the resampled texture has the expected size
        profile_array = nibabel.load(resampled_texturefile).get_data()
        profile_dim = profile_array.ndim
        profile_shape = profile_array.shape
        if profile_dim != 3:
            raise ValueError(
                "Expected profile texture array of dimension 3 not "
                "'{0}'".format(profile_dim))
        if (profile_shape[1] != 1) or (profile_shape[2] != 1):
            raise ValueError(
                "Expected profile texture array of shape (*, 1, 1) not "
                "'{0}'.".format(profile_shape))

        # Organize the resampled textures in a single file
        if sid in textures_map:
            raise ValueError("Subject '{0}' already treated. Check the intput "
                             "'regex'.".format(sid))
        textures_map[sid] = profile_array.ravel().astype(numpy.single)

    # Remove surfaces folder
    if not keep_individual_textures:
        shutil.rmtree(surfacesdir)

    # Save textures in CSV or in Numpy
    textures_files = []
    if save_mode in ["csv", "all"]:
        textures_file = os.path.join(fsoutdir, basename + ".csv")
        with open(textures_file, "wb") as open_file:
            csv_writer = csv.writer(open_file, delimiter=",")
            for sid in sorted(textures_map.keys()):
                row = [sid]
                row.extend(
                    ["{0:.4f}".format(elem) for elem in textures_map[sid]])
                csv_writer.writerow(row)
        textures_files.append(textures_file)
    if save_mode in ["numpy", "all"]:
        textures_file = os.path.join(fsoutdir, basename + ".npz")
        numpy.savez(textures_file, **textures_map)
        textures_files.append(textures_file)

    return textures_files


def mri_convert(fsdir, regex, output_directory=None, reslice=True,
                interpolation="interpolate",
                fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Export Freesurfer "*.mgz" image in Nifti format.

    Convert in native space: the destination image is resliced like the
    'rawavg.mgz' file if the reslice option is set. The converted file will
    then have a '.native' suffix.

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="regex" type="String" description="A regular expression
            used to locate the files to be converted from the 'fsdir'
            directory."/>
        <input name="output_directory" type="Directory" description="The
            conversion destination folder."/>
        <input name="reslice" type="Bool" description="If true reslice the
            input images like the raw image."/>
        <input name="interpolation" type="String" description="The
            interpolation method: interpolate|weighted|nearest|cubic."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="niftifiles" type="List_File" description="The converted
            nifti files."/>
    </unit>
    """
    # Check the interpolation method
    if interpolation not in ["interpolate", "weighted", "nearest", "cubic"]:
        raise ValueError(
            "'{0}' is not a valid interpolation method.".format(interpolation))

    # Get the images to convert from the regex
    inputs = glob.glob(os.path.join(fsdir, regex))
    if output_directory is not None:
        path = os.path.join(output_directory, "inputs.json")
        with open(path, "w") as open_file:
            json.dump(inputs, open_file, indent=4)

    # Convert each input file
    niftifiles = []
    for inputfile in inputs:

        # Create the output directory
        subject = inputfile.replace(fsdir, "")
        subject = subject.lstrip(os.sep).split(os.sep)[0]
        outdir = os.path.join(fsdir, subject, "convert")
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        # Create the FS command
        basename = os.path.basename(inputfile).replace(".mgz", "")
        cmd = ["mri_convert", "--resample_type", interpolation]
#               "--out_orientation", "RAS"]
        if reslice:
            reference_file = os.path.join(fsdir, subject, "mri", "rawavg.mgz")
            if not os.path.isfile(reference_file):
                raise ValueError("'{0}' does not exists, can't reslice image "
                                 "'{1}'.".format(reference_file, inputfile))
            cmd += ["--reslice_like", reference_file]
            basename = basename + ".native"
        converted_file = os.path.join(outdir, basename + ".nii.gz")
        niftifiles.append(converted_file)
        cmd += [inputfile, converted_file]

        # Execute the FS command
        recon = FSWrapper(cmd, shfile=fsconfig)
        recon()
        if recon.exitcode != 0:
            raise FreeSurferRuntimeError(
                recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr +
                recon.stdout)

    return niftifiles


def mri_vol2surf(hemi, volume_file, texture_file, ico_order, dat_file, fsdir,
                 sid,  surface_name="white",
                 fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Assigns values from a volume to each surface vertex.

    Wrapper around the FreeSurfer 'mri_vol2surf' command to create the
    described texture.

    Parameters
    ----------
    hemi: str (mandatory)
        hemisphere ('lh' or 'rh').
    volume_file: str (mandatory)
        input volume path.
    texture_file: str (mandatory)
        output texture file.
    ico_order: int (mandatory)
        icosahedron order in [0, 7] that will be used to generate the cortical
        surface texture at a specific tessalation (the corresponding cortical
        surface can be resampled using the
        'clindmri.segmentation.freesurfer.resample_cortical_surface' function).
    dat_file: str (mandatory)
        structural to FreeSurfer space affine transformation matrix '.dat'
        file as computed by 'tkregister2'.
    fsdir: str (mandatory)
        FreeSurfer subjects directory 'SUBJECTS_DIR'.
    sid: str (mandatory)
        FreeSurfer subject identifier.
    surface_name: str (optional, default 'white')
        The surface we  want to resample ('white' or 'pial').
    fsconfig: str (optional)
        The FreeSurfer '.sh' config file.
    """
    # Check input parameters
    if hemi not in ["lh", "rh"]:
        raise ValueError("'{0}' is not a valid hemisphere value which must be "
                         "in ['lh', 'rf']".format(hemi))
    if surface_name not in ["white", "pial"]:
        raise ValueError("'{0}' is not a valid surface value which must be in "
                         "['white', 'pial']".format(surface_name))
    if ico_order < 0 or ico_order > 7:
        raise ValueError("'Ico order '{0}' is not in 0-7 "
                         "range.".format(ico_order))

    # Construct the FS vol2surf command
    cmd = ["mri_vol2surf", "--src", volume_file, "--out", texture_file,
           "--srcreg", dat_file, "--hemi", hemi, "--trgsubject",
           "ico", "--icoorder", "{0}".format(ico_order), "--surf",
           surface_name, "--sd", fsdir, "--srcsubject", sid, "--noreshape",
           "--out_type", "mgz"]

    # Execute the FS command
    recon = FSWrapper(cmd, shfile=fsconfig)
    recon()
    if recon.exitcode != 0:
        raise FreeSurferRuntimeError(
            recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr + recon.stdout)


def mri_surf2surf(hemi, input_surface_file, output_surface_file, ico_order,
                  fsdir, sid,
                  fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Resample surface vertices.

    Wrapper around the FreeSurfer 'mri_surf2surf' command to create the
    described texture.

    Parameters
    ----------
    hemi: str (mandatory)
        hemisphere ('lh' or 'rh').
    input_surface_file: str (mandatory)
        input surface path.
    output_surface_file: str (mandatory)
        output '.mgz' surface path.
    ico_order: int (mandatory)
        icosahedron order in [0, 7] that will be used to generate the cortical
        surface texture at a specific tessalation (the corresponding cortical
        surface can be resampled using the
        'clindmri.segmentation.freesurfer.resample_cortical_surface' function).
    fsdir: str (mandatory)
        FreeSurfer subjects directory 'SUBJECTS_DIR'.
    sid: str (mandatory)
        FreeSurfer subject identifier.
    fsconfig: str (optional)
        The FreeSurfer '.sh' config file.
    """
    # Check input parameters
    if hemi not in ["lh", "rh"]:
        raise ValueError("'{0}' is not a valid hemisphere value which must be "
                         "in ['lh', 'rf']".format(hemi))
    if ico_order < 0 or ico_order > 7:
        raise ValueError("'Ico order '{0}' is not in 0-7 "
                         "range.".format(ico_order))

    # Set the output surface extension if necessary
    if not output_surface_file.endswith(".mgz"):
        output_surface_file += ".mgz"

    # Define FreeSurfer command
    cmd = ["mri_surf2surf", "--hemi", hemi, "--srcsurfval", input_surface_file,
           "--srcsubject", sid, "--trgsubject", "ico", "--trgicoorder",
           str(ico_order), "--trgsurfval", output_surface_file, "--sd", fsdir,
           "--trg_type", "mgz"]

    # Execute the FS command
    recon = FSWrapper(cmd, shfile=fsconfig)
    recon()
    if recon.exitcode != 0:
        raise FreeSurferRuntimeError(
            recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr + recon.stdout)


def resample_cortical_surface(
        fsdir, regex, output_directory=None, orders=[4, 5, 6, 7],
        surface_name="white",
        fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Resamples one cortical surface onto an icosahedron.

    Resample the white or pial FreeSurfer cotical surface using the
    'mri_surf2surf' command. Map also the associated annotation file.

    Can resample at different icosahedron order which specifies the size of the
    icosahedron according to the following table:
    Order  Number of Vertices
    0              12
    1              42
    2             162
    3             642
    4            2562
    5           10242
    6           40962
    7          163842

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="regex" type="String" description="A regular expression
            used to locate the surface files to be converted from the 'fsdir'
            directory."/>
        <input name="output_directory" type="Directory" description="The
            default resample destination folder."/>
        <input name="orders" type="List_Int" description="The icosahedron
            orders."/>
        <input name="surface_name" type="String" description="The surface we
            want to resample ('white' or 'pial')."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="resamplefiles" type="List_File" description="The
            resample surfaces."/>
        <output name="annotfiles" type="List_File" description="The
            resample annotations."/>
    </unit>
    """
    # Check input parameters
    if surface_name not in ["white", "pial"]:
        raise ValueError("'{0}' is not a valid surface value which must be in "
                         "['white', 'pial']".format(surface_name))
    norders = numpy.asarray(orders)
    if norders.min() < 0 or norders.max() > 7:
        raise ValueError("'At least one value in {0} is not in 0-7 "
                         "range.".format(orders))

    # Get all the subjects with the specified surface
    surfaces = glob.glob(os.path.join(fsdir, regex))
    if output_directory is not None:
        path = os.path.join(output_directory, "surfaces.json")
        with open(path, "w") as open_file:
            json.dump(surfaces, open_file, indent=4)

    # Go through all the subjects with the desired surface
    resamplefiles = []
    annotfiles = []
    for surf in surfaces:

        # Get some information based on the surface path
        subject_id = surf.split("/")[-3]
        hemi = os.path.basename(surf).split(".")[0]
        convertdir = os.path.join(fsdir, subject_id, "convert")
        if not os.path.isdir(convertdir):
            os.makedirs(convertdir)

        # Go through all specified orders
        for level in orders:

            # Construct the FS surface map command
            convertfile = os.path.join(convertdir, "{0}.{1}.{2}".format(
                hemi, surface_name, level))
            resamplefiles.append(convertfile)
            cmd = ["mri_surf2surf", "--sval-xyz", surface_name,
                   "--srcsubject", subject_id, "--trgsubject", "ico",
                   "--trgicoorder", str(level), "--tval", convertfile,
                   "--tval-xyz", "--hemi", hemi, "--sd", fsdir]

            # Execute the FS command
            recon = FSWrapper(cmd, shfile=fsconfig)
            recon()
            if recon.exitcode != 0:
                raise FreeSurferRuntimeError(
                    recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr +
                    recon.stdout)

            # Construct the FS label map command
            annotfile = os.path.join(convertdir, "{0}.aparc.annot.{1}".format(
                hemi, level))
            annotfiles.append(annotfile)
            if not os.path.isfile(annotfile):
                svalannot = os.path.join(fsdir, subject_id, "label",
                                         "{0}.aparc.annot".format(hemi))
                cmd = ["mri_surf2surf", "--srcsubject", subject_id,
                       "--trgsubject", "ico", "--trgicoorder", str(level),
                       "--hemi", hemi, "--sval-annot", svalannot,
                       "--tval", annotfile, "--sd", fsdir]

                # Execute the FS command
                recon = FSWrapper(cmd, shfile=fsconfig)
                recon()
                if recon.exitcode != 0:
                    raise FreeSurferRuntimeError(
                        recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr +
                        recon.stdout)

    # Remove duplicate annotation files
    annotfiles = list(set(annotfiles))

    return resamplefiles, annotfiles


def conformed_to_native_space(
        fsdir, regex, output_directory=None,
        fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Create a registration matrix between the conformed space (orig.mgz)
    and the native anatomical (rawavg.mgz).

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="regex" type="String" description="A regular expression
            used to locate the surface files to be converted from the 'fsdir'
            directory."/>
        <input name="output_directory" type="Directory" description="The
            default resample destination folder."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="trffiles" type="List_File"
            description="The conformed to native transformation files."/>
    </unit>
    """
    # Get all the subjects with a 'mri' directory
    mridirs = glob.glob(os.path.join(fsdir, regex))
    if output_directory is not None:
        path = os.path.join(output_directory, "mris.json")
        with open(path, "w") as open_file:
            json.dump(mridirs, open_file, indent=4)

    # Go through all the subjects with the desired folder
    trffiles = []
    for mdir in mridirs:

        # Get some information based on the folder path
        subject_id = mdir.rstrip("/").split("/")[-2]
        convertdir = os.path.join(fsdir, subject_id, "convert")
        if not os.path.isdir(convertdir):
            os.makedirs(convertdir)

        # Check that the two images of interest are present
        rawfile = os.path.join(mdir, "rawavg.mgz")
        origfile = os.path.join(mdir, "orig.mgz")
        if not (os.path.isfile(rawfile) and os.path.isfile(origfile)):
            raise ValueError("In folder '{0}' can't find file '{1}' or file "
                             "'{2}'.".format(mdir, rawfile, origfile))

        # Construct the FS command
        trffile = os.path.join(convertdir, "register.native.dat")
        trffiles.append(trffile)
        cmd = ["tkregister2", "--mov", rawfile, "--targ", origfile,
               "--reg", trffile, "--noedit", "--regheader"]

        # Execute the FS command
        recon = FSWrapper(cmd, shfile=fsconfig)
        recon()
        if recon.exitcode != 0:
            raise FreeSurferRuntimeError(
                recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr +
                recon.stdout)

    return trffiles


def surf_convert(fsdir, t1files, surffiles, output_directory=None,
                 rm_orig=False,
                 fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Export FreeSurfer surfaces to the native space.

    Note that all the vetices are given in the index coordinate system.
    The subjecy id in the t1 and surf files must appear in the -3 position:
        xxx/subject_id/convert/t1.nii.gz

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="output_directory" type="Directory" description="The
            conversion destination folder."/>
        <input name="t1files" type="List_File" description="The t1 nifti
            files."/>
        <input name="surffiles" type="List_File" description="The surface
            to be converted."/>
        <input name="rm_orig" type="Bool" description="If true remove
            the input surfaces."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="csurffiles" type="List_File" description="The converted
            surfaces in the native space."/>
    </unit>
    """
    # Create a t1 subject map
    t1map = {}
    for fname in t1files:
        subject_id = fname.split("/")[-3]
        if subject_id in t1map:
            raise ("Can't map two t1 for subject '{0}'.".format(subject_id))
        t1map[subject_id] = fname

    # Convert all the surfaces
    csurffiles = []
    for fname in surffiles:

        # Get the t1 reference image
        subject_id = fname.split("/")[-3]
        t1file = t1map[subject_id]
        t1_image = nibabel.load(t1file)

        # Compute the conformed space to the native anatomical deformation
        asegfile = os.path.join(fsdir, subject_id, "mri", "aseg.mgz")
        physical_to_index = numpy.linalg.inv(t1_image.get_affine())
        translation = tkregister_translation(asegfile, fsconfig)
        deformation = numpy.dot(physical_to_index, translation)

        # Load and warp the mesh
        # The mesh: a 2-uplet with vertex (x, y, z) coordinates and
        # mesh triangles
        mesh = freesurfer.read_geometry(fname)
        surf = TriSurface(vertices=apply_affine_on_mesh(mesh[0], deformation),
                          triangles=mesh[1])

        # Save the mesh in the native space
        outputfile = fname + ".native"
        surf.save(os.path.dirname(outputfile), os.path.basename(outputfile))
        csurffiles.append(outputfile)

        # Clean input surface if specified
        if rm_orig:
            os.remove(fname)

    return csurffiles


def qc(t1files, wmfiles, asegfiles, whitefiles, pialfiles, annotfiles,
       actor_ang=[0., 0., 0.], output_directory=None,
       fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Compute some quality check plots on the converted FrreSurfer
    outputs.

    The subjecy id in the input files must appear in the -3 position:
        xxx/subject_id/convert/t1.nii.gz

    Steps:

    * t1-images overlays
    * 3d surface segmentation snaps
    * t1-surfaces overlays

    actor_ang: float (optional, default 0)
        the actor rotation in the z direction.

    <unit>
        <input name="t1files" type="List_File" description="The
            t1 subject files."/>
        <input name="wmfiles" type="List_File" description="The
            white matter subject files."/>
        <input name="asegfiles" type="List_File" description="The
            subcortical segmentation subject files."/>
        <input name="output_directory" type="Directory" description="The
            conversion destination folder."/>
        <input name="whitefiles" type="List_File" description="The subject
            cortex surfaces."/>
        <input name="pialfiles" type="List_File" description="The subject pial
            surfaces."/>
        <input name="annotfiles" type="List_File" description="The pial/white
            surface annotations."/>
        <input name="actor_ang" type="List_Float" description="The actor x, y,
            z position (in degrees)."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="qcfiles" type="List_File" description="The quality check
            snaps."/>
    </unit>
    """
    import clindmri.plot.pvtk as pvtk
    from clindmri.plot.slicer import plot_image

    # Create a t1 subject map
    t1map = {}
    for fname in t1files:
        subject_id = fname.split("/")[-3]
        if subject_id in t1map:
            raise Exception("Can't map two t1 for subject '{0}'"
                            ".".format(subject_id))
        t1map[subject_id] = fname

    # Create the output list that will contain all the qc files
    qcfiles = []

    # Construct the t1-surfaces overlays and the 3d surface segmentation snaps
    ren = pvtk.ren()
    for name, files in [("white", whitefiles), ("pial", pialfiles)]:
        for fname in files:

            # Get the t1 reference image
            subject_id = fname.split("/")[-3]
            t1file = t1map[subject_id]
            t1_image = nibabel.load(t1file)

            # Get the qc output directory
            qcdir = os.path.join(os.path.dirname(fname), "qc")
            qcname = os.path.basename(fname)
            if not os.path.isdir(qcdir):
                os.makedirs(qcdir)

            # Get the triangular mesh
            basename = os.path.basename(fname).replace(
                name, "aparc.annot").replace(".native", "")
            annotfile = os.path.join(os.path.dirname(fname), basename)
            if annotfile not in annotfiles:
                raise ValueError(
                    "Annotation file '{0}' can't be found.".format(annotfile))
            surface = TriSurface.load(fname, annotfile=annotfile)

            # Construct the surfaces binarized volume
            binarizedfile = os.path.join(qcdir, qcname + ".nii.gz")
            overlay = numpy.zeros(t1_image.shape, dtype=numpy.uint)
            indices = numpy.round(surface.vertices).astype(int).T
            indices[0, numpy.where(indices[0] >= t1_image.shape[0])] = 0
            indices[1, numpy.where(indices[1] >= t1_image.shape[1])] = 0
            indices[2, numpy.where(indices[2] >= t1_image.shape[2])] = 0
            overlay[indices.tolist()] = 1
            overlay_image = nibabel.Nifti1Image(overlay, t1_image.get_affine())
            nibabel.save(overlay_image, binarizedfile)
            snap_file = os.path.join(qcdir, qcname + ".png")
            plot_image(t1file, overlay_file=binarizedfile, snap_file=snap_file,
                       name=qcname, overlay_cmap="cold_hot")
            qcfiles.append(snap_file)

            # Create a vtk surface actor of the cortex surface with the
            # associated labels
            ctab = [item["color"] for _, item in surface.metadata.items()]
            actor = pvtk.surface(
                surface.vertices, surface.triangles, surface.labels, ctab)
            actor.RotateX(actor_ang[0])
            actor.RotateY(actor_ang[1])
            actor.RotateZ(actor_ang[2])

            # Create a 3d surface segmentation snap
            pvtk.add(ren, actor)
            snaps = pvtk.record(ren, qcdir, qcname, n_frames=36,
                                az_ang=10, animate=True, delay=50)
            qcfiles.append(snaps[0])
            snaps = pvtk.record(ren, qcdir, qcname + ".3d", n_frames=1)
            qcfiles.append(snaps[0])
            pvtk.clear(ren)

    # Get the FreeSurfer lookup table
    fs_lut_names, fs_lut_colors = parse_fs_lut(os.path.join(
        os.path.dirname(fsconfig), "FreeSurferColorLUT.txt"))
    cmap = []
    nb_values = numpy.asarray(fs_lut_colors.keys()).max()
    cmap = numpy.zeros((nb_values, 4), dtype=numpy.single)
    for key, color in fs_lut_colors.items():
        if key > 0:
            cmap[key - 1, :3] = color
    cmap[:, 3] = 200.
    cmap /= 255.

    # Compute t1-images overlays
    for name, files in [("aseg", asegfiles), ("wm", wmfiles)]:
        for fname in files:

            # Get the t1 reference image
            subject_id = fname.split("/")[-3]
            t1file = t1map[subject_id]
            t1_image = nibabel.load(t1file)

            # Get the qc output directory
            qcdir = os.path.join(os.path.dirname(fname), "qc")
            if not os.path.isdir(qcdir):
                os.makedirs(qcdir)

            # Troncate the color map based on the label max
            array = nibabel.load(fname).get_data()
            order = sorted(set(array.flatten()))
            ccmap = cmap[order[1]: order[-1] + 1]

            # Overlay the current image with the t1 image
            qcname = "t1-{0}".format(name)
            snap_file = os.path.join(qcdir, qcname + ".png")
            plot_image(t1file, overlay_file=fname, snap_file=snap_file,
                       name=qcname, overlay_cmap=ccmap, cut_coords=(0, 0, 0))
            qcfiles.append(snap_file)

    return qcfiles


def parse_fs_lut(filename="FreeSurferColorLUT.txt"):
    """ Parse the FreeSurfer general lookup table.
    """
    fs_lut_names = {}
    fs_lut_colors = {}
    with open(filename) as open_file:
        for line in open_file:
            token = line.split()
            if len(token) >= 6:
                try:
                    fs_lut_names[int(token[0])] = token[1]
                    fs_lut_colors[int(token[0])] = (
                        int(token[2]), int(token[3]), int(token[4]))
                except ValueError:
                    continue
    return fs_lut_names, fs_lut_colors


def cortex(t1_file, fsdir, outdir, dest_file=None, prefix="cortex",
           generate_mask=True, generate_seeds=True):
    """ Compute a white matter mask and gyri labelization from the FreeSurfer
    'white' surface.

    Parameters
    ----------
    t1_file: str (mandatory)
        a file containing the t1 image used in FreeSurfer for the segmentation.
    fsdir: str( mandatory)
        the subject freesurfer segmentation directory.
    outdir: str (mandatory)
        the output directory.
    dest_file: str (optional, default None)
        a file containing an image where we want to project the segmentations:
        an affine transform is used to align this image to the t1 image.
    prefix: str (optional, default 'cortex')
        the output files prefix.
    generate_mask: bool (optional, default True)
        if True generate a white matter binary mask.
    generate_seeds: boll (optional, default False)
        if True create a 'seeds' directory containing all the gyri mask as
        idenpendent files.

    Returns
    -------
    mask_file: str
        the white matter mask image file.
    label_file: str
        the gyri label image file.
    seeds: list of str
        a list with the seed volumes.
    """
    # Create the output directory if necessary
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Load the dataset
    t1_image = nibabel.load(t1_file)
    t1_affine = t1_image.get_affine()

    # If a destination file is specified register it to the t1
    if dest_file is not None:

        # Load dataset
        dest_image = nibabel.load(dest_file)
        dest_affine = dest_image.get_affine()
        dest_shape = dest_image.get_shape()

        # In case of temporal serie extract the first volume
        if len(dest_shape) > 3:
            temporal_dest_file = dest_file
            dest_file = os.path.join(outdir, prefix + "_volume-0.nii.gz")
            extract_image(temporal_dest_file, index=0, out_file=dest_file)
            dest_shape = dest_shape[:3]

        # Register destination image to t1 image
        trf_file = os.path.join(outdir, prefix + "_dest_to_t1.trf")
        reg_file = os.path.join(outdir, prefix + "_dest_to_t1.nii.gz")
        flirt(dest_file, t1_file, omat=trf_file, out=reg_file, usesqform=False,
              cost="normmi", dof=6)
        voxel_dest_to_t1 = flirt2aff(trf_file, dest_file, t1_file)
        voxel_t1_to_dest = numpy.linalg.inv(voxel_dest_to_t1)

    # Otherwise use identity transformation
    else:
        trf_file = None
        reg_file = None
        dest_affine = t1_affine
        dest_shape = t1_image.get_shape()
        voxel_t1_to_dest = numpy.identity(4)

    # Load the FreeSurfer surface in the 'dest_file' voxel coordinates or
    # 't1_file' coordinates if not specified
    t1_physical_to_voxel = numpy.linalg.inv(t1_affine)
    seg = read_cortex_surface_segmentation(fsdir, t1_physical_to_voxel,
                                           voxel_t1_to_dest)

    # Create a mask of the white matter of both hemisphere
    if generate_mask:
        mask_array = seg["lh"].voxelize(dest_shape)
        mask_array += seg["rh"].voxelize(dest_shape)

    # Create a gyri label image of both hemisphere
    label_array = {}
    try:
        label_array["lh"], shift_lh = seg["lh"].labelize(dest_shape)
        label_array["rh"], shift_rh = seg["rh"].labelize(dest_shape, shift_lh)
    except:
        if reg_file is not None:
            raise FSLResultError("flirt")
        raise

    # Create the seeds
    seeds = []
    if generate_seeds:
        seedsdir = os.path.join(outdir, "gyri")
        if not os.path.isdir(seedsdir):
            os.mkdir(seedsdir)
        for hemi in ["lh", "rh"]:
            surf = seg[hemi]
            hemi_label_array = label_array[hemi]
            seed_array = numpy.zeros(hemi_label_array.shape,
                                     dtype=hemi_label_array.dtype)
            for index, item in surf.metadata.items():
                if index != 0:
                    if hemi == "rh":
                        index += shift_lh
                    seed_array[numpy.where(hemi_label_array == index)] = 1
                    seed_file = os.path.join(
                        seedsdir,
                        "{0}-{1}.nii.gz".format(hemi, item["region"]))
                    seed_image = nibabel.Nifti1Image(seed_array, dest_affine)
                    nibabel.save(seed_image, seed_file)
                    seed_array[...] = 0
                    seeds.append(seed_file)

    # Save the mask and label images
    mask_file = None
    if generate_mask:
        mask_file = os.path.join(outdir, prefix + "_mask.nii.gz")
        mask_image = nibabel.Nifti1Image(mask_array, dest_affine)
        nibabel.save(mask_image, mask_file)
    label_array = label_array["lh"] + label_array["rh"]
    label_file = os.path.join(outdir, prefix + "_gyri_labels.nii.gz")
    label_image = nibabel.Nifti1Image(label_array, dest_affine)
    nibabel.save(label_image, label_file)

    return mask_file, label_file, seeds, reg_file, trf_file


def population_statistic(fsdir, sid=None):
    """ Compute the mean std of each FreeSurfer individual scores. If sid
    is not None, return the subject values

    Parameters
    ----------
    fsdir: str (mandatory)
        the subject freesurfer segmentation directory.
    sid: str (optional, default None)
        if None compute the statistics on all subjects, otherwise return only
        the specified subject values.

    Returns
    -------
    popstats: dict
        the population statistics.
    """
    # Go through all statistics
    popstats = {
        "lh": {},
        "rh": {},
        "aseg": {}
    }
    stats = glob.glob(os.path.join(fsdir, "stats", "*.csv"))
    for fpath in stats:
        basename = os.path.basename(fpath).split(".")[0]
        if basename.startswith("aseg"):
            stype, _, sname = basename.split("_")
            hemi = "aseg"
            subject_header = "Measure:volume"
        elif basename.startswith("aparc"):
            stype, _, hemi, sname = basename.split("_")
            subject_header = "{0}.{1}.{2}".format(hemi, stype, sname)
        else:
            continue

        if sname not in popstats[hemi]:
            popstats[hemi][sname] = {}
        with open(fpath, "rb") as openfile:
            reader = csv.DictReader(openfile)
            for line in reader:
                subject = line.pop(subject_header)
                if sid is not None:
                    if subject != sid:
                        continue
                for key, value in line.items():
                    popstats[hemi][sname].setdefault(key, []).append(
                        float(value))
        for region_name, values in popstats[hemi][sname].items():
            mean = numpy.mean(values)
            std = numpy.std(values)
            popstats[hemi][sname][region_name] = {
                "values": values,
                "m": mean,
                "s": std
            }

    return popstats


class MutualExclusiveArg(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MandatoryArgMissing(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def mri_binarize(fsdir, sid, output_directory=None,
                 interpolation="interpolate",
                 region_ids=[], regions_label='-',
                 fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Generate binary mask from freesurfer segmentation.
        Regions can be defined either by their name (see list below) or their
        identification number

        wm
            Set match vals to 2, 41, 77, 251-255 (aseg for cerebral WM)

        ventricles
            Set match vals those for aseg ventricles+choroid (not 4th)

        wm+vcsf
            WM and ventricular CSF, including choroid (not 4th)

        subcort-gm
            Subcortical gray matter

        then the mask is converted in nifti format and reshaped to match the
        original file's dimensions.

    <unit>
        <input name="fsdir" type="Directory" description="The
            freesurfer working directory with all the subjects."/>
        <input name="sid" type="Str" description="The current subject
            identifier."/>
        <input name="output_directory" type="Directory" description="The
            freesurfer runtime folder."/>
        <input name="region_ids" type="List" content="Int" description="The
            identification numbers of the regions to segment"/>
        <input name="regions_label" type="Str" description="The region
            name, must be recognized by freesurfer"/>
        <input name="interpolation" type="Str" desc="Define the interpolation
            method, must be from 'interpolate', 'weighted', 'nearest', 'cubic'
            (default='interpolate')"/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="mask_file" type="File" description="The resulting mask"/>
    </unit>
    """
    # Create the fs output directory if necessary
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    # Check the interpolation method
    if interpolation not in ["interpolate", "weighted", "nearest", "cubic"]:
        raise ValueError(
            "'{0}' is not a valid interpolation method.".format(interpolation))

    # mutual exclusion check
    if len(region_ids) > 0 and regions_label != '-':
        raise MutualExclusiveArg("Both 'region_ids' and 'regions_label' are"
                                 " defined, only one of them should be")
    if not region_ids and not regions_label:
        raise MandatoryArgMissing("At least one of 'region_ids' or "
                                  "'regions_label' must be defined")

    fsfile = os.path.join(fsdir, sid, "mri", "aseg.mgz")
    reference_file = os.path.join(fsdir, sid, "mri", "rawavg.mgz")
    if not os.path.isfile(reference_file):
        raise ValueError("'{0}' does not exists, can't reslice image "
                         "'{1}'.".format(reference_file, fsfile))

    # Calls freesurfer: first generate binary mask then convert it in nifti
    # MASK generation
    # common part of the command
    if region_ids:
        mask_file = "freesurfer_{0}.mgz".format(
            "_".join([str(x) for x in region_ids]))
        out_mgz_path = os.path.join(output_directory, mask_file)
        cmd = ["mri_binarize", "--i", fsfile, "--o", out_mgz_path]
        for ids in region_ids:
            cmd.append("--match")
            cmd.append(str(ids))
    else:
        mask_file = "freesurfer_{0}.mgz".format(regions_label)
        out_mgz_path = os.path.join(output_directory, mask_file)
        regions = "--{0}".format(regions_label)
        cmd = ["mri_binarize", "--i", fsfile, "--o", out_mgz_path, regions]

    binarize = FSWrapper(cmd, shfile=fsconfig)
    binarize()
    if binarize.exitcode != 0:
        raise FreeSurferRuntimeError(
            binarize.cmd[0], " ".join(binarize.cmd[1:]),
            binarize.stderr + binarize.stdout)

    # CONVERSION
    mask_file = out_mgz_path.replace(".mgz", ".nii.gz")
    cmd = ["mri_convert", "--resample_type", interpolation,
           "--reslice_like", reference_file,
           out_mgz_path, mask_file]

    convert = FSWrapper(cmd, shfile=fsconfig)
    convert()
    if convert.exitcode != 0:
        raise FreeSurferRuntimeError(
            convert.cmd[0], " ".join(convert.cmd[1:]),
            convert.stderr + convert.stdout)

    return mask_file
