import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ne    = [15000, 12000, 10000, 8000, 6000, 4000, 2000, 1000, 800, 600, 400, 200, 150, 100, 80, 60, 40]
npack = [15000, 12000, 10000, 8000, 6000, 4000, 2000, 1000, 800, 600, 400, 200, 150, 100, 80, 60, 40]

# Arranged by [Npack] varies along rows vs [Ne] varies along columns

CPUtimeinmicrosec=np.array([[33452600, 21460100 , 14939100, 9521310, 5421300, 2366540, 588206, 155013, 100323, 56105, 25738, 6529, 3935, 1845, 1231, 730, 354],
                            [26697700, 17146000 , 11939600, 7606900, 4326480, 1885970, 469818, 124040, 79734 , 44516, 20586, 5267, 3172, 1499, 982 , 567, 271],
                            [22329900, 14330800 , 9986100 , 6365080, 3611690, 1570070, 393188, 103027, 66235 , 37348, 17255, 4416, 2656, 1255, 812 , 457, 223],
                            [17800900, 11419900 , 7951200 , 5070300, 2872440, 1247340, 313853, 82062 , 53041 , 29571, 13770, 3540, 2136, 1013, 630 , 356, 178],
                            [13397900, 8584590  , 5978800 , 3814100, 2154970, 934300 , 235123, 61394 , 39701 , 22359, 10367, 2663, 1602, 713 , 453 , 262, 134],
                            [8912950 , 5716930  , 3979840 , 2538630, 1426800, 619510 , 156699, 40888 , 26402 , 14927, 6904 , 1795, 1054, 462 , 300 , 176, 90 ],
                            [4464130 , 2858570  , 1992570 , 1264260, 710066 , 308243 , 78351 , 20512 , 13206 , 7471 , 3436 , 888 , 508 , 230 , 152 , 91 , 47 ],
                            [2256270 , 1443670  , 1002170 , 635995 , 356668 , 154084 , 39092 , 10334 , 6629  , 3737 , 1669 , 460 , 283 , 135 , 80  , 51 , 25 ],
                            [1826930 , 1169040  , 810763  , 515023 , 288621 , 123675 , 31283 , 8264  , 5318  , 2984 , 1338 , 362 , 225 , 105 , 62  , 41 , 20 ],
                            [1400900 , 891675   , 618717  , 392889 , 219912 , 93730  , 23614 , 6248  , 4005  , 2219 , 1004 , 273 , 171 , 81  , 48  , 32 , 21 ],
                            [967138  , 615835   , 427145  , 270916 , 152167 , 63628  , 16004 , 4198  , 2696  , 1468 , 667  , 184 , 115 , 53  , 31  , 20 , 13 ],
                            [436747  , 277852   , 192053  , 121360 , 67224  , 29911  , 7676  , 2025  , 1295  , 721  , 330  , 92  , 57  , 27  , 16  , 10 , 7  ],
                            [330391  , 209418   , 144601  , 91543  , 50798  , 22556  , 5789  , 1523  , 964   , 542  , 246  , 68  , 43  , 19  , 11  , 7  , 5  ],
                            [221596  , 140937   , 96977   , 61229  , 34287  , 15187  , 3903  , 1002  , 637   , 359  , 166  , 46  , 28  , 13  , 8   , 5  , 4  ],
                            [185910  , 113457   , 78549   , 49463  , 27910  , 12254  , 3130  , 801   , 524   , 295  , 135  , 36  , 23  , 11  , 7   , 5  , 4  ],
                            [145869  , 92889    , 63791   , 40857  , 22624  , 9323   , 2368  , 606   , 385   , 227  , 101  , 28  , 17  , 9   , 7   , 4  , 5  ],
                            [136729  , 87464    , 60902   , 38859  , 20088  , 6188   , 1593  , 413   , 265   , 148  , 69   , 20  , 12  , 6   , 5   , 5  , 3  ]])

GPUtimeinmicrosec=np.array([[24695700, 15497900 , 11037600 , 6847630 , 3802100 , 1690710, 430679 , 116998 , 74868.3, 41560.8, 25636.1, 8860.13, 6853.37, 6141.73, 6593.45, 6656.32, 6430.55],
                            [19739300, 12377000 , 8847730  , 5449720 , 3043340 , 1353470, 345558 , 93701.7, 59603  , 33292.4, 20437.5, 7130.75, 5417.71, 4834.35, 5254.04, 5292.34, 4969.91],
                            [16556700, 10320900 , 7419900  , 4539750 , 2536100 , 1128210, 287191 , 77734.8, 49675.7, 27690.9, 16738.7, 5805.98, 4419.44, 3967.82, 4246.19, 4219.33, 3997.41],
                            [13307000, 8241330  , 5984520  , 3619500 , 2027300 , 901204 , 229264 , 61854.6, 39663.6, 22085.8, 13201.9, 4545.26, 3484.41, 3055.78, 3276.61, 3175.24, 3073.57],
                            [10063500, 6371110  , 4524510  , 2701220 , 1519260 , 675413 , 171783 , 46334.2, 29793.9, 16597.3, 9673.55, 3343.13, 2520.79, 2215.52, 2296.38, 2265.42, 2200.34],
                            [6827290 , 4543710  , 3183990  , 1793980 , 1012170 , 449532 , 113980 , 30551  , 19552.7, 10981.5, 6094.95, 2054.35, 1544.22, 1264.55, 1302.24, 1294.21, 1258.09],
                            [3830150 , 2929100  , 2062480  , 894438  , 1305370 , 552295 , 56321.5, 14459.7, 9120.71, 5276.87, 2509.39, 782.137, 438.816, 395.507, 227.109, 149.653, 147.111],
                            [1792720 , 1009390  , 872933   , 447230  , 252285  , 111790 , 28244.5, 7291.07, 4578.23, 2661.6 , 1319.39, 419.724, 256.696, 224.521, 140.329, 106.096, 104.87 ],
                            [1268780 , 807965   , 562868   , 357908  , 201882  , 89445  , 22712.6, 5851.51, 3676.44, 2155.44, 1040.1 , 360.034, 210.588, 201.213, 131.162, 98.3345, 92.9058],
                            [950356  , 606888   , 422129   , 268620  , 151724  , 67261.9, 17038.5, 4372.83, 2786.74, 1640.8 , 801.583, 283.41 , 184.82 , 159.355, 111.382, 95.9915, 85.5452],
                            [634954  , 405293   , 281773   , 179547  , 101369  , 44935.5, 11391.9, 2979.01, 1900.16, 1115.11, 566.698, 215.781, 135.904, 123.93 , 89.595 , 78.4959, 86.4841],
                            [320017  , 204329   , 142334   , 90576.6 , 51156.7 , 22769.6, 5799.68, 1532.76, 1007.57, 598.727, 326.533, 136.987, 105.759, 91.7079, 85.07  , 70.394 , 67.0695],
                            [241698  , 154466   , 107677   , 68439.9 , 38668.8 , 17191.7, 4406.92, 1203.36, 768.261, 464.884, 258.133, 116.101, 85.7867, 84.1813, 73.3892, 67.6618, 68.6022],
                            [163011  , 104268   , 72503.1  , 46152.1 , 26033   , 11622.4, 3009.44, 817.235, 563.032, 332.478, 196.433, 97.2422, 87.4529, 73.578 , 64.2769, 63.1511, 64.8   ],
                            [131319  , 83953.5  , 58520.4  , 37201   , 21007.6 , 9366.57, 2416.13, 688.516, 444.837, 284.286, 175.055, 100.393, 73.2407, 69.6895, 61.2663, 64.3736, 63.0046],
                            [100249  , 64151.2  , 44720.3  , 28444.7 , 16139.9 , 7181.46, 1871.14, 541.056, 369.032, 233.393, 137.224, 80.4498, 71.8693, 71.6241, 66.4725, 63.3763, 59.9208],
                            [69013.8 , 44133.1  , 30774.1  , 19575.6 , 11139.9 , 4946.32, 1300.79, 397.474, 267.981, 182.707, 117    , 74.2817, 68.616 , 67.9105, 61.0163, 67.7712, 59.6102]])


CPUgflopspersec=np.zeros(np.shape(CPUtimeinmicrosec))
GPUgflopspersec=np.zeros(np.shape(GPUtimeinmicrosec))

for idx, i in enumerate(npack):
    print "jdx,j=",idx,i
    for jdx, j in enumerate(ne):
        print "jdx,j=",jdx,j
        CPUgflopspersec[idx, jdx] = (i * j*j * 1e-6) / (CPUtimeinmicrosec[idx, jdx]*1.0e-3)
        GPUgflopspersec[idx, jdx] = (i * j*j * 1e-6) / (GPUtimeinmicrosec[idx, jdx]*1.0e-3)

print "GPU=",GPUgflopspersec

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)
    ax.set_ylabel(r'N$_{pack}$', fontsize=16)
    ax.set_title(r'N$_{e}$', fontsize=16)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, fontsize=12, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontsize=10)
    ax.set_yticklabels(row_labels, fontsize=10)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


fig, ax = plt.subplots()
im, cbar = heatmap(CPUgflopspersec, ne, npack, ax=ax,
                   cmap="gist_rainbow", cbarlabel="GFLOP/s")
fig.tight_layout()
plt.savefig('FFM_CPU.png', dpi=300)

fig, ax = plt.subplots()
im, cbar = heatmap(GPUgflopspersec, ne, npack, ax=ax,
                   cmap="gist_rainbow", cbarlabel="GFLOP/s")
fig.tight_layout()
plt.savefig('FFM_GPU.png', dpi=300)


