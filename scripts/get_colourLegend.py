"""quick script to put a colour legend on an annotated microgaph png"""

# call with: python3 get_colourLegend.py ../data/ img9.flat.spots.png 20 0 FALSE ibm
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
    maxdiam = sys.argv[3]  # float/nm
    mindiam = sys.argv[4]  # float/nm, is this always 0?
    noscale = sys.argv[5].lower() == "true"  # bool, if true put scale bar in pixels
    colourscaleuse = sys.argv[6].lower()

    outputname = str("colourbarred_" + ImageName)

    imageIn = read_image(WorkDirectory, ImageName)
    print(imageIn.shape)
    imgdims = imageIn.shape  # generalise later

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

    ax.set_yticks([])
    ax.set_xticks([])

    if colourscaleuse == "rainbow":
        nucmap = mpl.colors.LinearSegmentedColormap.from_list(
            "",
            [
                [0.5, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 0.5, 0.0],  # matplot lib needs RGB vals as floats
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.75, 0.75],
                [0.0, 0.0, 1.0],
                [0.75, 0.0, 0.75],
                [1.0, 1.0, 1.0],
            ],
        )
    elif colourscaleuse == "ibm":
        nucmap = mpl.colors.LinearSegmentedColormap.from_list(
            "",
            [
                [0.39215686, 0.56078431, 1.0],
                [0.43137255, 0.46666667, 0.97254902],
                [
                    0.47058824,
                    0.36862745,
                    0.94117647,
                ],  # matplot lib needs RGB vals as floats
                [0.66666667, 0.25882353, 0.72156863],
                [0.8627451, 0.14901961, 0.49803922],
                [0.92941176, 0.26666667, 0.24705882],
                [0.99607843, 0.38039216, 0.0],
                [1.0, 0.5372549, 0.0],
                [1.0, 0.69019608, 0.0],
            ],
        )
    else:
        print("Arg 6, colourscale, is unkown deafulting to IBM")
        nucmap = mpl.colors.LinearSegmentedColormap.from_list(
            "",
            [
                [0.39215686, 0.56078431, 1.0],
                [0.43137255, 0.46666667, 0.97254902],
                [
                    0.47058824,
                    0.36862745,
                    0.94117647,
                ],  # matplot lib needs RGB vals as floats
                [0.66666667, 0.25882353, 0.72156863],
                [0.8627451, 0.14901961, 0.49803922],
                [0.92941176, 0.26666667, 0.24705882],
                [0.99607843, 0.38039216, 0.0],
                [1.0, 0.5372549, 0.0],
                [1.0, 0.69019608, 0.0],
            ],
        )
    if noscale:
        fig.colorbar(
            mpl.cm.ScalarMappable(
                norm=mpl.colors.Normalize(mindiam, maxdiam), cmap=nucmap
            ),
            ax=ax,
            orientation="vertical",
            label="Loop Diameter/px",
        )
    else:
        fig.colorbar(
            mpl.cm.ScalarMappable(
                norm=mpl.colors.Normalize(mindiam, maxdiam), cmap=nucmap
            ),
            ax=ax,
            orientation="vertical",
            label="Loop Diameter/nm",
        )

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
