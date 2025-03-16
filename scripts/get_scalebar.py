"""quick script to put a scale bar on a microgaph png"""

# call with: python get_scalebar.py ../data/ img9.png 0.51 F
import sys
import os
import skimage
from skimage import io
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib as mpl
from matplotlib.pyplot import savefig


def read_image(fpath, fname):
    """Function to read image.

    Performs a check on that the file format is supported, a check for
    file existance is perfomred by another fuction (CheckImageExists)

    Parameters
    ----------
    Fpath : str
        The path of the file to read
    fname : str
        The name of the file including extension

    Returns
    -------
    testpic1 : ndarray
            MxNxC array, dtype type dependent on file format

    Raises
    ------
    NameError
        If the file format is not .jpg, .png or .tiff
    """
    # Split the extension from the path and normalise it to lowercase.
    ext = os.path.splitext(fname)[-1].lower()
    if ext in [".jpg", ".png", ".tiff"]:
        fname = os.path.join(fpath, fname)
        testpic1 = skimage.io.imread(fname)
        return testpic1
    else:
        raise NameError(str(ext) + " is not a supported filetype")


def main():
    """Main body function"""
    WorkDirectory = str(sys.argv[1])
    ImageName = str(sys.argv[2])
    pixellength = float(sys.argv[3])  # float/nm
    negativeIN = sys.argv[
        4
    ]  # boolean, true -> image is brightfield so use black text, false -> image is darkfield so use white text

    outputname = str("scalebarred_" + ImageName)
    if negativeIN == "T":
        negative = True
    elif negativeIN == "F":
        negative = False
    else:
        negative = True
        print("WARN: argument 4 invalide input (T or F expected), defaulting to T")

    if negative:
        print("DBG negative ==true, black text")
        fontcolour = "black"  # could allow white box with black text if mixed later.
    else:
        print("DBG negative ==false, white text")
        fontcolour = "white"

    imageIn = read_image(WorkDirectory, ImageName)
    print(imageIn.shape)
    imgdims = imageIn.shape  # generalise later

    image_length = imgdims[1] * pixellength

    fontprops = fm.FontProperties(size=18)
    fig, ax = plt.subplots()
    print(imgdims[0])
    print(imgdims[1])
    ax.imshow(
        imageIn,
        extent=[0, (imgdims[1] - 1), 0, (imgdims[0] - 1)],
        cmap="gray",
        interpolation="none",
    )

    roughlen = image_length / 5  # in pixels
    roughlen = roughlen * pixellength  # in nm
    print(f"rough len is {roughlen}")
    base = 50
    targetLength = base * round(roughlen / base)
    print(f"targetLength is {targetLength}")
    targetLengthStr = str(str(targetLength) + " nm")
    targetLength = targetLength / pixellength  # convert back to pixels for drawing

    scalebar = AnchoredSizeBar(
        ax.transData,
        targetLength,
        targetLengthStr,
        "lower right",  #'upper left', 'upper center', 'upper right', 'center left',
        pad=0.1,  #'center', 'center right', 'lower left', 'lower center', 'lower right'
        color=fontcolour,
        frameon=False,
        size_vertical=5,
        fontproperties=fontprops,
    )

    ax.add_artist(scalebar)
    ax.set_yticks([])
    ax.set_xticks([])
    # -------save image-----
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    dirout = str("../data/" + outputname)
    print(f"DBG out put dir is: {dirout}")
    savefig(dirout, transparent=True, dpi=600)  # edgecolor='none', facecolor='none'


if __name__ == "__main__":
    main()
