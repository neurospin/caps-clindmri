##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os

# Clindmri import
from clindmri.extensions.fsl import FSLWrapper
from clindmri.extensions.fsl.exceptions import FSLRuntimeError


def flirt(in_file, ref_file, omat=None, out=None, init=None, cost="corratio",
          usesqform=None, displayinit=None, anglerep="euler", bins=256,
          interp="trilinear", dof=12, applyxfm=None, verbose=0,
          shfile="/etc/fsl/5.0/fsl.sh"):
    """ Wraps command flirt.

    Usage: flirt [options] -in <inputvol> -ref <refvol> -out <outputvol>
           flirt [options] -in <inputvol> -ref <refvol> -omat <outputmatrix>
           flirt [options] -in <inputvol> -ref <refvol> -applyxfm -init
               <matrix> -out <outputvol>

    Available options are:
        -in  <inputvol>
            (no default)
        -ref <refvol>
            (no default)
        -init <matrix-filname>
            (input 4x4 affine matrix)
        -omat <matrix-filename>
            (output in 4x4 ascii format)
        -out, -o <outputvol>
            (default is none)
        -datatype {char,short,int,float,double}
            (force output data type)
        -cost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}
            (default is corratio)
        -searchcost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}
            (default is corratio)
        -usesqform
            (initialise using appropriate sform or qform)
        -displayinit
            (display initial matrix)
        -anglerep {quaternion,euler}
            (default is euler)
        -interp {trilinear,nearestneighbour,sinc,spline}
            (final interpolation: def - trilinear)
        -sincwidth <full-width in voxels>
            (default is 7)
        -sincwindow {rectangular,hanning,blackman}
        -bins <number of histogram bins>
            (default is 256)
        -dof  <number of transform dofs>
            (default is 12)
        -noresample
            (do not change input sampling)
        -forcescaling
            (force rescaling even for low-res images)
        -minsampling <vox_dim>
            (set minimum voxel dimension for sampling (in mm))
        -applyxfm
            (applies transform (no optimisation) - requires -init)
        -applyisoxfm <scale>
            (as applyxfm but forces isotropic resampling)
        -paddingsize <number of voxels>
            (for applyxfm: interpolates outside image by size)
        -searchrx <min_angle> <max_angle>
            (angles in degrees: default is -90 90)
        -searchry <min_angle> <max_angle>
            (angles in degrees: default is -90 90)
        -searchrz <min_angle> <max_angle>
            (angles in degrees: default is -90 90)
        -nosearch
            (sets all angular search ranges to 0 0)
        -coarsesearch <delta_angle>
            (angle in degrees: default is 60)
        -finesearch <delta_angle>
            (angle in degrees: default is 18)
        -schedule <schedule-file>
            (replaces default schedule)
        -refweight <volume>
            (use weights for reference volume)
        -inweight <volume>
            (use weights for input volume)
        -wmseg <volume>
            (white matter segmentation volume needed by BBR cost function)
        -wmcoords <text matrix>
            (white matter boundary coordinates for BBR cost function)
        -wmnorms <text matrix>
            (white matter boundary normals for BBR cost function)
        -fieldmap <volume>
            (fieldmap image in rads/s - must be already registered to the
            reference image)
        -fieldmapmask <volume>
            (mask for fieldmap image)
        -pedir <index>
            (phase encode direction of EPI - 1/2/3=x/y/z & -1/-2/-3=-x/-y/-z)
        -echospacing <value>
            (value of EPI echo spacing - units of seconds)
        -bbrtype <value>
            (type of bbr cost function: signed [default], global_abs,
            local_abs)
        -bbrslope <value>
            (value of bbr slope)
        -setbackground <value>
            (use specified background value for points outside FOV)
        -noclamp
            (do not use intensity clamping)
        -noresampblur
            (do not use blurring on downsampling)
        -2D
            (use 2D rigid body mode - ignores dof)
        -verbose <num>
            (0 is least and default)
        -v
            (same as -verbose 1)
        -i
            (pauses at each stage: default is off)
        -version
            (prints version number)
        -help
    """
    # Set default parameters
    dirname = os.path.dirname(in_file)
    basename = os.path.basename(in_file).split(".")[0]
    if out is None:
        out = os.path.join(dirname, "flirt_out_{0}.nii.gz".format(basename))
    if omat is None:
        omat = os.path.join(dirname, "flirt_omat_{0}".format(basename))

    # Call flirt
    fslprocess = FSLWrapper("flirt", shfile=shfile)
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]),
                              fslprocess.stderr)

    return out, omat
