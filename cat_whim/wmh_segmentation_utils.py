# Define imports
import numpy as np
import nibabel as nib

from pathlib import Path
from cat_whim.config import INTERIM_DATA_DIR, BIDS_DATA_DIR, INCL_EXCL_FLW_DIR


TMP_DIR = Path("/tmp")
WMH_DIR = BIDS_DATA_DIR / "derivatives" / "wmh_segmentations"

if not Path.exists(WMH_DIR):
    Path.mkdir(WMH_DIR)

def clear_tmp_dir(TMP_DIR):
    """
    Clears the temporary directory of your system from niftis. This is done because ANTS does not this automatically, or,
    at least, I wasn't able to do this automatically with ANTS..
    """
    # Ants is making a lot of temporary files that do not get automatically removed...
    for file in TMP_DIR.iterdir():
        # Check if the file starts with "tmp"
        if file.is_file() and file.name.endswith(".nii.gz"):
            # Delete the file
            file.unlink()

def create_box_mni6():
    # Create a 3D numpy array filled with ones
    array = np.ones((182, 218, 182), dtype=float)

    # Define the indices to set to zero. This was done empirically to cover the area near the olfactory bulbs
    # where most of the wrong segmentation was observed.
    x_indices = slice(75, 107) 
    y_indices = slice(130, 180) 
    z_indices = slice(30, 60) 

    # Set the values to zero within the specified range
    array[x_indices, y_indices, z_indices] = 0.
    return array

def from_nibabel(nib_image):
    """
    Convert a nibabel image to an ANTsImage.
    Note that AntsPy0.5.3 removed ants.to_nibabel, so this is the same code used before, taken from:
    https://antspyx.readthedocs.io/en/latest/_modules/ants/utils/convert_nibabel.html#from_nibabel
    """
    import os
    from tempfile import mkstemp
    from ants import ants_image_io as iio2
    fd, tmpfile = mkstemp(suffix=".nii.gz")
    nib_image.to_filename(tmpfile)
    new_img = iio2.image_read(tmpfile)
    os.close(fd)
    os.remove(tmpfile)
    return new_img


def to_nibabel(image):
    """
    Convert an ANTsImage to a Nibabel image. 
    Note that AntsPy0.5.3 removed ants.to_nibabel, so this is the same code used before, taken from:
    https://antspyx.readthedocs.io/en/latest/_modules/ants/utils/convert_nibabel.html#to_nibabel
    """
    import nibabel as nib
    import os
    from tempfile import mkstemp

    fd, tmpfile = mkstemp(suffix=".nii.gz")
    image.to_filename(tmpfile)
    new_img = nib.load(tmpfile)
    os.close(fd)
    # os.remove(tmpfile) ## Don't remove tmpfile as nibabel lazy loads the data.
    return new_img

def postproc_wmh_mask_mni6(wmh_mask_mni6):
    filterbox = create_box_mni6()
    wmh_mask_mni6_data = wmh_mask_mni6.get_fdata()
    # Multiply the box filter array by the wmh mask data
    result = filterbox * wmh_mask_mni6_data
    # Get the affine matrix from the original image
    affine = wmh_mask_mni6.affine
    # Create the new nibabel image with the result array and the same affine matrix
    result_img = nib.Nifti1Image(result, affine)
    return result_img