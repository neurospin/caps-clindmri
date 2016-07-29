##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import glob

# Clindmri import
from clindmri.extensions.fsl import FSLWrapper
from clindmri.extensions.fsl.exceptions import FSLRuntimeError


def probtrackx2(samples, seed, mask, dir, stop=None, avoid=None,
                out="fdt_paths", nsamples=5000, nsteps=2000, cthr=0.2,
                loopcheck=None, waypoints=None, onewaycondition=None,
                usef=None, simple=None, seedref=None, steplength=0.5,
                fibthresh=0.01, distthresh=0.0, sampvox=0.0, network=False,
                omatrix1=False, omatrix2=False, omatrix3=False,
                target2=None, target3=None,
                shfile="/etc/fsl/5.0/fsl.sh"):
    """ Wraps command probtrackx2.

    Single voxel
    ------------

    [1] Connectivity from a single seed point.

    probtrackx2(samples="/.../fsl.bedpostX/merged",
                mask="/.../fsl.bedpostX/nodif_brain_mask",
                seed="$PATH/tracto/seedvox_coordinates.txt",
                simple=True,
                loopcheck=True,
                dir="$PATH",
                out="SingleVoxel_paths")

    [2] Tracking in a standard / no-diffusion space.

    probtrackx2(samples="/.../fsl.bedpostX/merged",
                mask="/.../fsl.bedpostX/nodif_brain_mask",
                seed="$PATH/tracto/seedvox_coordinates.txt",
                seeref="/.../fsl.bedpostX/nodif_brain_mask",
                simple=True,
                loopcheck=True,
                dir="$PATH",
                out="SingleVoxel_paths")

    Single mask
    -----------

    probtrackx2(seed="/.../lh-precentral.nii.gz",
                loopcheck=True,
                onewaycondition=True,
                samples="/.../fsl.bedpostX/merged",
                mask="/.../fsl.bedpostX/nodif_brain_mask",
                dir="$PATH")

    Multiple masks
    --------------
    probtrackx2(network=True,
                seed="$PATH/masks.txt",
                loopcheck=True,
                onewaycondition=True,
                samples="/.../fsl.bedpostX/merged",
                mask="/.../fsl.bedpostX/nodif_brain_mask",
                dir="$PATH")

    Usage:
    probtrackx2 -s <basename> -m <maskname> -x <seedfile> -o <output>
                --targetmasks=<textfile>
    probtrackx2 --help

    Compulsory arguments (You MUST set one or more of):
        -s,--samples Basename for samples files - e.g. 'merged'
        -m,--mask    Bet binary mask file in diffusion space
        -x,--seed    Seed volume or list (ascii text file) of volumes and/or
                     surfaces

    Optional arguments (You may optionally specify one or more of):
        -o,--out	   Output file (default='fdt_paths')
        --dir       Directory to put the final volumes in - code makes
                    this directory - default='logdir'
        --forcedir  Use the actual directory name given - i.e. don't add +
                    to make a new directory
        --simple    Track from a list of voxels (seed must be a ASCII list of
                    coordinates)
        --network   Activate network mode - only keep paths going through at
                    least one of the other seed masks
        --opd       Output path distribution
        --pd        Correct path distribution for the length of the pathways
        --fopd      Other mask for binning tract distribution
        --os2t      Outputs seeds to target images. One per voxel in the seed.
                    There can be quite a lot of these files.
        --s2tastext Output seed-to-target counts as a text file (default in
                    simple mode)
        --targetmasks File containing a list of target masks - for
                      seeds_to_targets classification
        --waypoints Waypoint mask or ascii list of waypoint masks - only keep
                    paths going through ALL the masks
        --waycond   Waypoint condition. Either 'AND' (default) or 'OR'
        --wayorder  Reject streamlines that do not hit waypoints in given
                    order. Only valid if waycond=AND
        --onewaycondition Apply waypoint conditions to each half tract
                          separately
        --avoid     Reject pathways passing through locations given by
                    this mask
        --stop      Stop tracking at locations given by this mask file
        --omatrix1  Output matrix1 - SeedToSeed Connectivity
        --distthresh1	Discards samples (in matrix1) shorter than this
                      threshold (in mm - default=0)
        --omatrix2  Output matrix2 - SeedToLowResMask
        --target2   Low resolution binary brain mask for storing connectivity
                    distribution in matrix2 mode
        --omatrix3  Output matrix3 (NxN connectivity matrix)
        --target3   Mask used for NxN connectivity matrix (or Nxn if
                    lrtarget3 is set)
        --lrtarget3 Column-space mask used for Nxn connectivity matrix
        --distthresh3 Discards samples (in matrix3) shorter than this
                      threshold (in mm - default=0)
        --omatrix4  Output matrix4 - DtiMaskToSeed (special Oxford Sparse
                    Format)
        --colmask4  Mask for columns of matrix4 (default=seed mask)
        --target4   Brain mask in DTI space
        --xfm       Transform taking seed space to DTI space (either FLIRT
                    matrix or FNIRT warpfield) - default is identity
        --invxfm	   Transform taking DTI space to seed space (compulsory when
                    using a warpfield for seeds_to_dti)
        --seedref   Reference vol to define seed space in simple mode -
                    diffusion space assumed if absent
        --meshspace Mesh reference space - either 'caret' (default) or
                    'freesurfer' or 'first' or 'vox'

        -P,--nsamples   Number of samples - default=5000
        -S,--nsteps     Number of steps per sample - default=2000
        --steplength    Steplength in mm - default=0.5

        --distthresh    Discards samples shorter than this threshold (in mm -
                        default=0)
        -c,--cthr       Curvature threshold - default=0.2
        --fibthresh     Volume fraction before subsidary fibre orientations are
                        considered - default=0.01
        -l,--loopcheck  Perform loopchecks on paths - slower, but allows lower
                        curvature threshold
        -f,--usef       Use anisotropy to constrain tracking
        --modeuler  Use modified euler streamlining

        --sampvox       Sample random points within x mm sphere seed voxels
                        (e.g. --sampvox=5). Default=0
        --randfib   Default 0. Set to 1 to randomly sample initial fibers
                    (with f > fibthresh).
                    Set to 2 to sample in proportion fibers
                    (with f>fibthresh) to f.
                    Set to 3 to sample ALL populations at random
                    (even if f<fibthresh)
        --fibst     Force a starting fibre for tracking - default=1, i.e.
                    first fibre orientation. Only works if randfib==0
        --rseed     Random seed

    Parameters
    ----------
    shfile: str (optional, default NeuroSpin path)
        the path to the FSL 'fsl.sh' configuration file.

    Returns
    -------
    proba_files: list of str
        a list of files containing probabilistic fiber maps.
    network_file: str
        a voxel-by-target connection matrix.
    """
    # Call probtrackx
    fslprocess = FSLWrapper("probtrackx2 --opd --forcedir", shfile=shfile,
                            optional="ALL")
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]),
                              fslprocess.stderr)

    # Get the outputs
    proba_files = glob.glob(os.path.join(dir, out + "*"))
    network_file = None
    if network is not None:
        network_file = os.path.join(dir, "fdt_network_matrix")

    return proba_files, network_file


def bedpostx(input, n=3, w=1, b=1000, j=1250, s=25, model=2, g=None, c=None,
             rician=None, shfile="/etc/fsl/5.0/fsl.sh", cpus=""):
    """ Wraps command bedpostx.

    Usage: bedpostx <input> [options]

    expects to find bvals and bvecs in input
    expects to find data and nodif_brain_mask in input
    expects to find grad_dev in input, if -g is set

    options (old syntax)

    -n (number of fibers per voxel, default 3)
    -w (ARD weight, more weight means less secondary fibers per voxel,
       default 1)
    -b (burnin period: number of iterations before starting the sampling.
        These might be increased if the data are noisy, and the MCMC needs more
        iterations to converge, default 1000)
    -j (number of jumps to be made by MCMC, default 1250)
    -s (sample every, default 25)
    -model (deconvolution model.
           1: single-shell, with sticks,
           2: multi-shell, with sticks with a range of diffusivities (default),
           3: multi-shell, with zeppelins)
    -g (consider gradient nonlinearities,
       instructs bedpostx to use the grad_dev.nii.gz file from the data folder
       and produce voxel-specific bvals and bvecs, default off)
    --rician (use a Rician noise modeling to replace the default Gaussian
               noise assumption)
    -c do not use CUDA capable hardware/queue (if found)

    ALTERNATIVELY: you can pass on xfibres options onto directly bedpostx
    For example:  bedpostx <subject directory> --noard --cnonlinear
    Type 'xfibres --help' for a list of available options
    Default options will be bedpostx default (see above), and not xfibres
    default.

    Parameters
    ----------
    shfile: str (optional, default NeuroSpin path)
        the path to the FSL 'fsl.sh' configuration file.

    Returns
    -------
    outdir: str
        The bedpostx output directory
    merged_th<i>samples - 4D volume
        Samples from the distribution on theta
    merged_ph<i>samples - 4D volume
        Samples from the distribution on phi: theta and phi together represent
        the principal diffusion direction in spherical polar co-ordinates
    merged_f<i>samples - 4D volume
        Samples from the distribution on anisotropic volume fraction.
    mean_th<i>samples - 3D Volume
        Mean of distribution on theta
    mean_ph<i>samples - 3D Volume
        Mean of distribution on phi
    mean_f<i>samples - 3D Volume
        Mean of distribution on f anisotropy. Note that in each voxel, fibers
        are ordered according to a decreasing mean f-value
    mean_dsamples - 3D Volume
        Mean of distribution on diffusivity d
    mean_S0samples - 3D Volume
        Mean of distribution on T2w baseline signal intensity S0
    dyads<i>
        Mean of PDD distribution in vector form. Note that this file can be
        loaded into fslview for easy viewing of diffusion directions
    dyads_dispersion - 3D Volume
        Uncertainty on the estimated fiber orientation. Characterizes how wide
        the orientation distribution is around the respective PDD.
    nodif_brain_mask
        binary mask created from nodif_brain - copied from input directory
    bvecs
        contain a 3x1 vector for each gradient, indicating the gradient
        direction - copied from input directory
    bvals
        contain a scalar value for each applied gradient, corresponding to the
        respective bvalue - copied from input directory
    """
    # Call bedpostx
    fslprocess = FSLWrapper("bedpostx", shfile=shfile, cpus=cpus)
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]),
                              fslprocess.stderr)

    # Format outputs
    outdir = input + ".bedpostX"
    merged_th = glob.glob(os.path.join(outdir, "merged_th*"))
    merged_ph = glob.glob(os.path.join(outdir, "merged_ph*"))
    merged_f = glob.glob(os.path.join(outdir, "merged_f*"))
    mean_th = glob.glob(os.path.join(outdir, "mean_th*"))
    mean_ph = glob.glob(os.path.join(outdir, "mean_ph*"))
    mean_f = glob.glob(os.path.join(outdir, "mean_f*"))
    mean_d = os.path.join(outdir, "mean_d*")
    mean_S0 = os.path.join(outdir, "mean_S0*")
    dyads = glob.glob(os.path.join(outdir, "dyads*"))

    return (outdir, merged_th, merged_ph, merged_f, mean_th, mean_ph, mean_f,
            mean_d, mean_S0, dyads)


def bedpostx_datacheck(input, shfile="/etc/fsl/5.0/fsl.sh"):
    """ Wraps bedpostx_datacheck

    Usage: bedpostx_datacheck input

    Parameters
    ----------
    input: str (mandatory)
        the folder to check.
    shfile: str (optional, default NeuroSpin path)
        the path to the FSL 'fsl.sh' configuration file.

    Returns
    -------
    is_valid: bool
        True if all the data are present in the input directory
    """
    # Call bedpostx_datacheck
    fslprocess = FSLWrapper("bedpostx_datacheck", shfile=shfile)
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]),
                              fslprocess.stderr)

    # Parse outputs
    is_valid = (
        fslprocess.stderr == "" and "does not exist" not in fslprocess.stdout)

    return is_valid
