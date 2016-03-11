# -*- coding: utf-8 -*-

import subprocess
import os
import shutil


def freesurfer_snaps_wm(fsdir, sid, axis=["C", "A", "S"],
                        slice_interval=[0, 255],
                        output_directory=None,
                        fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh",
                        logfile=None):
    """
    plot snaps from freesurfer segmentation, using freesurfer tkmedit

    """

    # sanity checks
    for item in axis:
        if item not in ["C", "A", "S"]:
            raise Exception("Axis '{}' is unknown, use 'C' for Coronal, 'A' "
                            "for Axial and 'S' for sagittal")

    axis_orient = {"C": 0,
                   "A": 1,
                   "S": 2}

    # set env viariable
    os.environ['SUBJECTS_DIR'] = fsdir

    try:
        if not os.path.isfile(os.path.join(fsdir, sid, "mri", "nu.mgz")):
            raise Exception("{}: nu.mgz file missing".format(sid))

        # create output dir
        if os.path.isdir(os.path.join(output_directory, sid)):
            shutil.rmtree(os.path.join(output_directory, sid))
        os.makedirs(os.path.join(output_directory, sid))

        # generate the tcl script
        script_path = os.path.join(output_directory, "wm_pial.tcl")
        with open(script_path, "w") as script:
            script.write("LoadMainSurface 0 lh.white\n")
            script.write("LoadMainSurface 1 rh.white\n")

            for axe in axis:
                script.write("SetCursor 0 128 128 128\n")
                script.write("SetOrientation {}\n".format(axis_orient[axe]))
                script.write("SetCursor 0 0 0 0\n")
                script.write("for { set slice %s } { $slice <= %s } "
                             "{ incr slice 1 } {\n" % (slice_interval[0],
                                                       slice_interval[1]))
                script.write("\tSetSlice $slice\n")
                script.write("\tRedrawScreen\n")
                script.write("\tSaveRGB {}" .format(
                    os.path.join(output_directory, sid,
                                 "snapshot-wm-{}-$slice.rgb\n".format(axe))))
                script.write("}\n\n")

            script.write("UnloadSurface 0\n")
            script.write("UnloadSurface 1\n")
            script.write("QuitMedit\n")

        subprocess.check_call(["tkmedit",
                               sid,
                               "nu.mgz",
                               "-tcl",
                               script_path])
        os.remove(script_path)

        # converting rgb files
        for snap in os.listdir(os.path.join(output_directory, sid)):
            subprocess.check_call([
                "convert",
                os.path.join(output_directory, sid, snap),
                os.path.join(output_directory, sid,
                             "{}.png".format(snap.split(".")[0]))])
            os.remove(os.path.join(output_directory, sid, snap))
    except:
        if logfile is not None:
            if os.path.isfile(logfile):
                arg = "a"
            else:
                arg = "w"
            with open(logfile, arg) as _file:
                _file.write("{}\n".format(sid))
        else:
            pass
