

# System import
import os
import re
import numpy
import subprocess
import nibabel
import scipy.sparse
from nibabel import freesurfer
from nibabel.gifti import read, write, GiftiDataArray, GiftiImage
import matplotlib.pyplot as plt
import clindmri.plot.pvtk as pvtk

# Dipy import
from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from dipy.tracking.eudx import EuDX
from dipy.reconst import peaks, shm
from dipy.tracking import utils


def mesh_connectivity(mesh, mask=None):
    """ Generate a connectivity matrix from a mesh and apply a region masking.
    """
    # Unpack mesh information
    vertex, triangles = mesh
    vertex_x, vertex_y, vertex_z = vertex.T

    mask = None
    # Build mesh edges
    edges = []
    if mask == None:
        mask = numpy.arange(len(vertex_x))
        new_idx = numpy.arange(len(vertex_x)) 
    else:
        new_idx = numpy.cumsum(mask) - 1
    for triangle in triangles:
        a, b, c = triangle
        a_, b_, c_ = new_idx[triangle]
        if (mask[a] and mask[b]):
            edges.append([a_, b_])
            edges.append([b_, a_])
        if (mask[c] and mask[b]):
            edges.append([b_, c_])
            edges.append([c_, b_])
        if (mask[a] and mask[c]):
            edges.append([a_, c_])
            edges.append([c_, a_])
                    
    edges = numpy.array(edges)
    print edges
    connectivity = scipy.sparse.coo_matrix((numpy.ones(edges.shape[0]),
                                   (edges.T[0], edges.T[1])))
    return connectivity


def physical_to_index(vertex, physical_to_index_matrix, translation=None,
                      round=True):
    """ Transform a physical coordinate to an index coordinate.

    If 'round' is 'True', the index will be rounded and cast to the nearest 
    integer value, otherwise it will be kept as a floating-point value.
    """
    if translation is None:
        translation = numpy.eye(4)
    transformation = numpy.dot(physical_to_index_matrix, translation)
    ones = numpy.ones((vertex.shape[0], 1))
    homogenous_vertex = numpy.concatenate((vertex, ones), axis=1)
    indices = numpy.dot(transformation, homogenous_vertex.T).T
    if round:
        indices = indices.round().astype(int)
    return indices[:, :3]


def label_image(mesh, labels, index_to_physical_matrix, translation_vector,
                label_array, shift=0):
    """ Create an image from the mesh and associated labels.
    """
    physical_to_index_matrix = numpy.linalg.inv(index_to_physical_matrix)
    translation = numpy.eye(4)
    translation[:3, 3] = translation_vector
    indices = physical_to_index(mesh[0], physical_to_index_matrix, translation)
    nb_of_labels = 0
    for label in set(labels):
        if label != -1:
            nb_of_labels += 1
            label_indices = indices[numpy.where(labels == label)]
            label_array[label_indices.T.tolist()] = label + shift
    return nb_of_labels


def get_ras(mgz_aseg_file, tkregister=False):
    """ Access the image coordiante space.
  
    Freesurfer use a special origin for the Right-Anterior-Superior
    (anatomical coordinates) space: to get the the scanner Vox2RAS 
    'mri_info --vox2ras aseg.mgz' vs the FreeSurfer/tkregister Vox2RAS
    'mri_info --vox2ras-trk aseg.mgz'.
    """
    command = ["mri_info", "--vox2ras", mgz_aseg_file]
    if tkregister:
        command[1] = "--vox2ras-tkr"
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    returncode = process.returncode
    stdout = ",".join([line.strip() for line in stdout.splitlines()])
    stdout = re.sub(r"  *", ",", stdout)
    ras = numpy.fromstring(stdout, dtype=float, sep=",").reshape(4, 4)

    return ras


def label_to_diffusion(t1file, labelsfile, dfile, trffile, anatlabelfile):
    """ Register the t1 image to the diffusion image and warp the labels
    to the diffusion space.
    """
    # Affine registration
    command = ["flirt", "-in", t1file, "-ref", dfile, "-dof", "12", "-omat",
           trffile, "-cost", "normmi"]
    subprocess.check_call(command)

    # Warp labels
    command = ["flirt", "-in", labelsfile, "-ref", dfile, "-out", anatlabelfile,
               "-init", trffile, "-applyxfm", "-interp", "nearestneighbour"]
    subprocess.check_call(command)
    


def supervised_anatomical_connectivity(dfile, bvecfile, bvalfile, labelsfile,
                                       order=4, density=1, step=0.5):
    """ 
    """
    # Read diffusion sequence
    bvals, bvecs = read_bvals_bvecs(bvalfile, bvecfile)
    gtab = gradient_table(bvals, bvecs)
    ddata = nibabel.load(dfile).get_data()

    # Read the labels
    labels = nibabel.load(labelsfile).get_data()

    # Estimate ODF model
    white_matter = labels == 0
    csamodel = shm.CsaOdfModel(gtab, order)
    csapeaks = peaks.peaks_from_model(
        model=csamodel, data=ddata, sphere=peaks.default_sphere,
        relative_peak_threshold=.8, min_separation_angle=45,
        mask=white_matter)

    # Compute deterministic tractography
    seeds = utils.seeds_from_mask(white_matter, density=density)
    streamline_generator = EuDX(
        csapeaks.peak_values, csapeaks.peak_indices,
        odf_vertices=peaks.default_sphere.vertices, a_low=.05, step_sz=step,
        seeds=seeds)
    affine = streamline_generator.affine
    streamlines = list(streamline_generator)

    # Compute the connectivity matrix
    connectivity, grouping = utils.connectivity_matrix(
        streamlines, labels, affine=affine, return_mapping=True,
        mapping_as_streamlines=True)

    # Remove the background
    connectivity = connectivity[1:, 1:]

    return connectivity

    
wmfile = "/volatile/imagen/dmritest/000000022453/fs/mri/wm.seg.mgz"    
dfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.nii.gz"
bvecfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.bvec"
bvalfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.bval"
fsdir = "/volatile/imagen/dmritest/000000022453/fs"
asegfile = "/volatile/imagen/dmritest/000000022453/fs/mri/aseg.mgz"
t1file = "/volatile/imagen/dmritest/000000022453/ADNI_MPRAGE/000000022453s012a1001.nii.gz"
outdir = "/volatile/imagen/dmritest/000000022453/processed"
mask_label_name = "postcentral"
meshdir = os.path.join(fsdir, "surf")
labeldir = os.path.join(fsdir, "label")
labelsfile = os.path.join(outdir, "labels.nii.gz")
trffile = os.path.join(outdir, "anat_to_diff.trf")
anatlabelfile = os.path.join(outdir, "anat_labels.nii.gz")


# Qunatize the cortex surface segmentation
# Load the t1 image
t1im = nibabel.load(t1file)
label_array = numpy.zeros(t1im.get_shape(), dtype=numpy.int16)
shift = 0
for side in ["lh", "rh"]:
    # Load the mesh: a 2-uplet with vertex (x, y, z) coordinates and
    # mesh triangles
    meshfile = os.path.join(meshdir, "{0}.white".format(side))
    mesh = freesurfer.read_geometry(meshfile)

    # Get annotation id at each vertex (if a vertex does not belong
    # to any label and orig_ids=False, its id will be set to -1) and
    # the names of the labels
    annotfile = os.path.join(labeldir, "{0}.aparc.annot".format(side))
    labels, ctab, names = freesurfer.read_annot(annotfile, orig_ids=False)

    # Plot test
    ren = pvtk.ren()
    actor = pvtk.surface(mesh[0], mesh[1], labels, ctab)
    pvtk.add(ren, actor)
    pvtk.show(ren)

    # Get the translation betwwen the ras/ras-trk spaces
    ras = get_ras(asegfile, tkregister=False)
    rastrk = get_ras(asegfile, tkregister=True)
    translation = (ras - rastrk)[:3, 3]

    # Create the label mask
    shift = label_image(mesh, labels, t1im.get_affine(), translation,
                        label_array, shift)

    # Create a connectivity matrix describing the mesh structure
    #if mask_label_name in names:
    #    mask = labels == names.index("postcentral")
    #else:
    #    mask = None
    #connectivity = mesh_connectivity(mesh, mask=mask)
    #print connectivity

# Save the image label
labelim = nibabel.Nifti1Image(label_array, t1im.get_affine())
nibabel.save(labelim, labelsfile)

# Get segmentation
#command = ["mri_convert", wmfile, 

# Register the anatomical and dissufion sequance
label_to_diffusion(t1file, labelsfile, dfile, trffile, anatlabelfile)

# Compute the connectivity matrix
aconnectivity = supervised_anatomical_connectivity(
    dfile, bvecfile, bvalfile, anatlabelfile, order=4, density=1, step=0.5)
plt.imshow(numpy.log1p(aconnectivity), interpolation="nearest")
plt.show()

