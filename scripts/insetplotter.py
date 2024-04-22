import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.colors
import numpy as np
from operator import add
import os
import argparse
import enum
import sys
import json
import re

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


sentinalStr = "SeNtInAl"
tmpFile = "plotlines.tmp"
maxColumns = 7

# Plotter cheat sheet
#
# XY
#   SeNtInAl, tag_type, user_defined, tag_name, x, y
#
#
# BAR | GROUPED_BAR | BARBOX | BOX
#   SeNtInAl, tag_type, user_defined, tag_name, x_label, y
#
#       tag_names - plotting multiple will create segmented bars with
#                   a legend. Legend can be overridden with --custom-legend-labels
#       x_label - string used to denote X
#       y       - float value
#
# THREEDBAR
#   SeNtInAl, tag_type, z, tag_name, x_label, y_label
#       tag_names - plotting multiple will create segmented bars with
#                   a legend. Legend can be overridden with --custom-legend-labels
#       x_label - string used to denote X
#       y_label - string used to denote Y
#       z       - float



class TAG_TYPE(enum.Enum):
    UNDEFINED = 0
    XY = 1
    HISTOGRAM = 2
    CDF = 3
    BAR = 4
    BOX = 5
    BARBOX = 6
    GROUPED_BAR = 7
    THREEDBAR = 8
    THREEDSURF = 9
    THREEDPROJECTX = 10
    THREEDPROJECTY = 11
    THREEDUNROLL = 12

#class dracula():
#   CYAN = '#8be9fd'
#   GREEN = '#50fa7b'
#   ORANGE = '#ffb86c'
#   PINK = '#ff79c6'
#   PURPLE = '#bd93f9'
#   RED = '#ff5555'
#   YELLOW = '#f1fa8c'
#   COLORS = [CYAN, ORANGE, GREEN, PINK, PURPLE, YELLOW, RED]

class dracula():
   ZERO = '#D291AC'
   ONE = '#779E44'
   TWO = '#F6B7AA'
   THREE = '#367E51'
   FOUR = '#5F1333'
   FIVE = '#426615'
   SIX = '#AD5B4A'
   #COLORS = [ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, ZERO, ONE, TWO, THREE, FOUR]
   COLORS = [plt.cm.Pastel1(i) for i in range(9)]

#cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", [dracula.PINK, dracula.PURPLE])
#cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", [dracula.CYAN, dracula.ORANGE])

cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", [dracula.ZERO, dracula.TWO])
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", [dracula.FOUR, dracula.SIX])

#cms = [cm.inferno, cm.viridis]
#cms = [cm.inferno, cm.cool]
#cms = [cm.copper, cm.summer]
cms = [cmap1, cmap2]


class tag:
    def __init__(self):
        self.tagType = TAG_TYPE.UNDEFINED
        self.name = ""



def getTags(csvFilePaths, tagWhitelist):
    tags = []

    printlogwarning = True;
    with open(tmpFile, "w") as aggregateFile:
        for p in csvFilePaths:
            with open(p, "r") as csvIn:
                lines = csvIn.readlines()

            for line in lines:
                if line[:len(sentinalStr)] == sentinalStr:
                    tokens = line.split(',')

                    if len(tokens) < 6 or len(tokens) > 8:
                        if printlogwarning:
                            print("WARNING - Malformed log line")
                            print(line)
                            printlogwarning = False
                        continue

                    # strip trailing \n, add right amount of commas, add \n back
                    line = line.rstrip() + (','*(maxColumns - len(tokens))) + '\n'
                    aggregateFile.write(line)

                    if tokens[3] == "error":
                        print("WARNING - OVERFLOW ERROR IN LOG FILE")
                        print(line)
                        continue

                    whitelisted = tagWhitelist == None
                    if not whitelisted:
                        for whitelistedTag in tagWhitelist:
                            whitelisted = whitelisted or re.match('^'+whitelistedTag+'$', tokens[3])
                    if not whitelisted:
                        continue

                    newTag = True
                    for t in tags:
                        if t.name == tokens[3]: # must be exact match to be old tag
                            newTag = False
                            break

                    if newTag and whitelisted:
                        print("found new tag: " + tokens[3])
                        t = tag()
                        if tokens[1] == "xy":
                            t.tagType = TAG_TYPE.XY
                        if tokens[1] == "box":
                            t.tagType = TAG_TYPE.BOX
                        elif tokens[1] == "histogram":
                            t.tagType = TAG_TYPE.HISTOGRAM
                        elif tokens[1] == "cdf":
                            t.tagType = TAG_TYPE.CDF
                        elif tokens[1] == "bar":
                            t.tagType = TAG_TYPE.BAR
                        elif tokens[1] == "grouped_bar":
                            t.tagType = TAG_TYPE.GROUPED_BAR
                        elif tokens[1] == "barbox":
                            t.tagType = TAG_TYPE.BARBOX
                        elif tokens[1] == "3dbar":
                            t.tagType = TAG_TYPE.THREEDBAR
                        elif tokens[1] == "3dsurf":
                            t.tagType = TAG_TYPE.THREEDSURF
                        elif tokens[1] == "3dprojectx":
                            t.tagType = TAG_TYPE.THREEDPROJECTX
                        elif tokens[1] == "3dprojecty":
                            t.tagType = TAG_TYPE.THREEDPROJECTY
                        elif tokens[1] == "3dunroll":
                            t.tagType = TAG_TYPE.THREEDUNROLL
                        t.name = tokens[3]
                        tags.append(t)

    # make sure tags in same order as tagWhitelist
    if tagWhitelist != None:
        newTags = []
        for tname in tagWhitelist:
            for t in tags:
                if re.match('^'+tname+'$', t.name):
                    newTags.append(t)
                    continue
        return newTags
    else:
        return tags



def plot(csvFilePaths, options, tags):
    data = np.genfromtxt(tmpFile, dtype=None, invalid_raise = False,
            delimiter=',', encoding=None, filling_values="0",
            names="sentinal, type, func, tag, x, y, z")

    tagIndex = 0;
    tagLabels = []
    lastPlotBottom = []
    legendArtistProxyShapes=[]
    skipLegend = options.hide_legend

    for t in tags:
        tagLabels.append(t.name)

    plt.figure(0, figsize=(options.fig_w, options.fig_h), dpi=800)

    for t in tags:
        print("plotting tag: " + t.name)

        if t.tagType == TAG_TYPE.XY:
            exes = []
            whys = []
            for row in data:
                if row[3] == t.name and (options.xfilter == None or row[4] in options.xfilter):
                    exes.append(float(row[4]));
                    whys.append(float(row[5]));
            # decorate, sort, then undecorate
            exes, whys = (np.array(t) for t in zip(*sorted(zip(exes, whys))))
            # average duplicate values
            dedup_exes = np.unique(exes)
            dedup_whys = np.empty(dedup_exes.shape)
            for i, x in enumerate(dedup_exes):
                dedup_whys[i] = np.mean(whys[exes == x]) # / 1e9

            if options.title == "Network IO per Iteration":
              dedup_whys /= 1e9

            if t.name == 'emp_float_lm_time_vs_points_per_loc_itr_latency':
                edgecolor=dracula.COLORS[tagIndex+options.color_offset]
                import colorsys
                h, l, s = colorsys.rgb_to_hls(*matplotlib.colors.to_rgb(edgecolor))
                edgehighlight = colorsys.hls_to_rgb(h, min(1, l * 0.5), s = s)

                plt.plot(dedup_exes, dedup_whys,
                    label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                    color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None),
                    linestyle = '--', lw=2, path_effects=[pe.Stroke(linewidth=3, foreground=edgehighlight), pe.Normal()])
            else:
                plt.plot(dedup_exes, dedup_whys,
                    label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                    color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))

        elif t.tagType == TAG_TYPE.HISTOGRAM:
            print("plotting histogram tag: " + t.name)
            print("TODO!")

        elif t.tagType == TAG_TYPE.BAR or t.tagType == TAG_TYPE.GROUPED_BAR or t.tagType == TAG_TYPE.BARBOX or t.tagType == TAG_TYPE.BOX:
            labels = [] # vector of x-axis label strings
            dataByLabel = [] # vector of all y values, indexed by label
            offsetDataByLabel = [] # dataByLable, stacked
            whys = [] # average y value indexed by label
            width = 0.5  # width of a group of bars for grouped bar graph

            # collect labels for all tags
            for row in data:
                if len(row) < 6:
                    print("bar|grouped_bar|box|barbox data bad format")
                    exit()
                for tagLabel in tagLabels:
                    if row[3] == tagLabel:
                        if row[4] not in labels and (options.xfilter == None or row[4] in options.xfilter):
                            labels.append(row[4])
                            dataByLabel.append([])
                            offsetDataByLabel.append([])

            if len(lastPlotBottom) == 0:
                lastPlotBottom = [0]*len(labels)

            # collect data for this tag
            for row in data:
                if row[3] == t.name and (options.xfilter == None or row[4] in options.xfilter):
                    try:
                        i = labels.index(row[4])
                        dataByLabel[i].append(float(row[5]))
                        offsetDataByLabel[i].append(float(row[5]) + lastPlotBottom[i])
                    except:
                        continue
            for row in dataByLabel:
                if len(row) != 0:
                    whys.append(np.average(row))
                else:
                    whys.append(0)

            if t.tagType == TAG_TYPE.BOX:
                # no label on box plots
                plt.boxplot(dataByLabel, positions=range(len(labels)))#, showmeans=True)

            if t.tagType == TAG_TYPE.BAR:
                if options.horizontal:
                    plt.barh(range(len(whys)), whys, bottom=lastPlotBottom,
                        label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))
                else:
                    plt.bar(range(len(whys)), whys, bottom=lastPlotBottom,
                        label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))
                if len(tags) > 1 and not skipLegend:
                    plt.legend() # bar always gets legend, must be after plotting
                lastPlotBottom = list( map(add, lastPlotBottom, whys) )# for stacking using average

            if t.tagType == TAG_TYPE.GROUPED_BAR:
                offset=len(tags)/2*width/5
                if options.horizontal:
                    rects = plt.barh(np.arange(len(whys)) - offset + (width*tagIndex/len(tags)), whys, width/len(tags),
                        label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None),
                        edgecolor='black', linewidth=1, hatch=options.hatching)
                else:
                    # special case for snail paper
                    if t.name == "EMP_MUL":
                        rects = plt.bar(np.arange(len(whys)) - offset + (width*tagIndex/len(tags)), whys, width/len(tags),
                            label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                            color="#2ca02c",
                            edgecolor='black', linewidth=1, hatch=options.hatching)
                    else:
                        rects = plt.bar(np.arange(len(whys)) - offset + (width*tagIndex/len(tags)), whys, width/len(tags),
                            label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                            color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None),
                            edgecolor='black', linewidth=1, hatch=options.hatching)
                if len(tags) > 1 and not skipLegend:
                    plt.legend() # bar always gets legend, must be after plotting
                lastPlotBottom = list( map(add, lastPlotBottom, whys) )# for stacking using average

                # label rectangles
                if options.bar_label_dict and t.name in options.bar_label_dict:
                    if len(options.bar_label_dict[t.name]) != len(whys):
                        print('Error - wrong number of bar labels - %d (should be %d)\n\n\n' % (len(options.bar_label_dict[t.name]), len(whys)) )
                        return
                    rectnum=0
                    for rect in rects:
                        height = rect.get_height()
                        plt.text(rect.get_x() + rect.get_width()/2., 1.05*height, options.bar_label_dict[t.name][rectnum],
                                ha='center', va='bottom',
                                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
                        rectnum=rectnum+1


            if t.tagType == TAG_TYPE.BARBOX:
                # stack boxplot data - shift upwards by last avg
                meanpointprops = dict(marker='o', markeredgecolor='black', markerfacecolor='black')
                r = plt.boxplot(offsetDataByLabel, positions=range(len(labels)), showmeans=True, meanprops=meanpointprops)

                # need to parse returned means because they dont include outliers
                whysExcludingOutliers = []
                for i in range(len(labels)):
                    whysExcludingOutliers.append(r['means'][i].get_ydata()[0] - lastPlotBottom[i])

                plt.bar(range(len(whysExcludingOutliers)), whysExcludingOutliers, width, bottom=lastPlotBottom,
                    label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                    color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None),
                    edgecolor='black', linewidth=1, hatch=options.hatching)
                lastPlotBottom = list( map(add, lastPlotBottom, whys) )# for stacking using average
                if len(tags) > 1 and not skipLegend:
                    plt.legend()

            if options.print_raw:
                print(labels)
                np.set_printoptions(threshold=sys.maxsize)
                print(whys)

            if options.rotate_x_labels:
                plt.xticks(range(len(labels)), labels, rotation=45)
            else:
                if options.horizontal:
                    plt.yticks(range(len(labels)), labels)
                else:
                    plt.xticks(range(len(labels)), labels)


        elif t.tagType == TAG_TYPE.THREEDBAR or t.tagType == TAG_TYPE.THREEDSURF or t.tagType == TAG_TYPE.THREEDPROJECTX or t.tagType == TAG_TYPE.THREEDPROJECTY or t.tagType == TAG_TYPE.THREEDUNROLL:
            if t.tagType == TAG_TYPE.THREEDBAR or t.tagType == TAG_TYPE.THREEDSURF:
                ax = plt.gca(projection = '3d')

            #plt.rcParams['legend.fontsize'] = 10
            xlabels = [] # vector of x-axis label strings
            ylabels = [] # vector of y-axis label strings

            for row in data:
                if len(row) < 7:
                    print("3d bar data bad format")
                    exit()
                if row[3] == t.name:
                    if (row[4] not in xlabels and (options.xfilter == None or row[4] in options.xfilter)):
                        xlabels.append(row[4])
                    if (row[5] not in ylabels):
                        ylabels.append(row[5])

            if len(lastPlotBottom) == 0:
                lastPlotBottom = [0]*len(xlabels)*len(ylabels)
            elif len(lastPlotBottom) != len(xlabels)*len(ylabels):
                print("ERROR - tags have different number of labels")
                print(len(lastPlotBottom))
                print(xlabels)
                print(ylabels)
                exit()

            zDataByLabel = [[[] for x in range(len(ylabels))] for y in range(len(xlabels))] # vector of all z values, indexed by x/y label
            zees = [[0 for x in range(len(ylabels))] for y in range(len(xlabels))] # average z value indexed by x/y labels

            for row in data:
                if row[3] == t.name and (options.xfilter == None or row[4] in options.xfilter):
                    xi = xlabels.index(row[4])
                    yi = ylabels.index(row[5])
                    zDataByLabel[xi][yi].append(float(row[6]))

            for i in range(len(xlabels)):
                for j in range(len(ylabels)):
                    zRow = zDataByLabel[i][j]
                    if len(zRow) != 0:
                        zees[i][j] = np.average(zRow)

            x = []
            y = []
            dz = []
            for i in range(len(xlabels)):
                for j in range(len(ylabels)):
                    x.append(i);
                    y.append(j);
                    dz.append(zees[i][j]);
            dx = np.ones(len(x))
            dy = np.ones(len(y))

            if options.print_raw:
                print(xlabels)
                print(ylabels)
                np.set_printoptions(threshold=sys.maxsize)
                print(np.resize(dz, [len(xlabels), len(ylabels)]))

            ### Note - this should work:
            ###     surf = ax.bar3d(x, y, lastPlotBottom, dx, dy, dz, label=t.name)
            ### But it looks terrible!
            ### Plotting the 0 height bars breaks the rendering.
            ### Now lets remove them carefully so we don't break lastPlotBottom

            filteredx = []
            filteredy = []
            filteredlpb = []
            filtereddz = []
            for i in range(len(xlabels)):
                for j in range(len(ylabels)):
                    if zees[i][j] != 0:
                        filteredx.append(i)
                        filteredy.append(j)
                        filtereddz.append(zees[i][j])
                        filteredlpb.append(lastPlotBottom[i*len(ylabels)+j])
            filtereddx = np.ones(len(filtereddz))
            filtereddy = np.ones(len(filtereddz))

            if t.tagType == TAG_TYPE.THREEDBAR:
                surf = ax.bar3d(filteredx, filteredy, filteredlpb, filtereddx, filtereddy, filtereddz,
                        label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))
                #surf = ax.bar3d(x, y, lastPlotBottom, dx, dy, dz, label=t.name)

                # Default legend looks terrible
                #surf._facecolors2d=surf._facecolors3d
                #surf._edgecolors2d=surf._edgecolors3d
                #plt.legend()
                # Use this instead: https://stackoverflow.com/questions/5803015/how-to-create-a-legend-for-3d-bar-in-matplotlib
                legendArtistProxyShapes.append(plt.Rectangle((0,0), 1, 1,
                    fc=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else plt.get_cmap("tab10")(tagIndex))))
                ax.set_zlabel(options.zlabel)
                if options.azimuth:
                    ax.azim = options.azimuth
                if options.elevation:
                    ax.elev = options.elevation

                plt.xticks(range(len(xlabels)), xlabels)#, rotation=45)
                plt.yticks(range(len(ylabels)), ylabels)#, rotation=45)

            elif t.tagType == TAG_TYPE.THREEDSURF:
                surf = ax.plot_trisurf(x, y, dz, cmap=cms[tagIndex+options.color_offset], linewidth=0,
                        antialiased=False, alpha=(tagIndex+1)/len(tags))
                clb = plt.colorbar(surf, shrink=0.5, aspect=14,
                        pad=(.16 if (tagIndex == len(tags)-1) else .04))
                clb.ax.set_title(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name)
                # set_xlabel('hi') to put at bottom
                ax.set_zlabel(options.zlabel)
                plt.xticks(range(len(xlabels)), xlabels)#, rotation=45)
                plt.yticks(range(len(ylabels)), ylabels)#, rotation=45)
                plt.tight_layout()# dont move
                skipLegend = True

            elif t.tagType == TAG_TYPE.THREEDPROJECTX:
                if options.projection == None:
                    print("need projection command line option")
                    return
                else:
                    sliced = []
                    slicedLPB = []
                    for i in range(len(xlabels)):
                        sliced.append(zees[i][options.projection])
                        slicedLPB.append(lastPlotBottom[i * len(ylabels) + options.projection])

                    # Plot a line
                    plt.plot([float(i) for i in xlabels], sliced,
                        label=options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name,
                        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))

                    # Stacked Bars
                    #plt.bar([i for i in xlabels], sliced, bottom=slicedLPB,
                    #        label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                    #        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))

                    plt.xticks([float(i) for i in xlabels], xlabels)#, rotation=45)

            elif t.tagType == TAG_TYPE.THREEDPROJECTY:
                if options.projection == None:
                    print("need projection command line option")
                    return
                else:
                    # Plot a line
                    plt.plot([float(i) for i in ylabels], zees[options.projection],
                        label=options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name,
                        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))

                    # Stacked Bars
                    #plt.bar([i for i in xlabels], zees[options.projection], bottom=lastPlotBottom[options.projection*len(ylabels) : options.projection*len(ylabels)+len(ylabels)], # THIS MAY NOT BE RIGHT
                    #        label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                    #        color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))

                    plt.xticks([float(i) for i in ylabels], ylabels)#, rotation=45)

            elif t.tagType == TAG_TYPE.THREEDUNROLL:
                width = .3 # width of bar
                numX = len(xlabels)


                # insert spaces
                spacedLPB = lastPlotBottom.copy()
                spacedDz = dz.copy()
                for i in range(len(ylabels)):
                    if i != 0 and i != len(ylabels):
                        spacedDz.insert((i*len(ylabels)) + i - 1, 0)
                        spacedLPB.insert((i * len(ylabels)) + i - 1, 0)

                plt.bar(range(len(spacedDz)), spacedDz, bottom=spacedLPB,
                    label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                    color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))
                if len(tags) > 1 and not skipLegend:
                    plt.legend() # bar always gets legend, must be after plotting

                xylabels = []
                for i in range(len(xlabels)):
                    if i != 0:
                        xylabels.append('')
                    for j in range(len(ylabels)):
                        xylabels.append(str(xlabels[i]) + ' - ' + str(int(ylabels[j])))
                plt.xticks(range(len(xylabels)), xylabels, rotation=45)

            lastPlotBottom = list( map(add, lastPlotBottom, dz) )# for stacking using average


        tagIndex += 1;

    if len(legendArtistProxyShapes) > 1 and not skipLegend: # 3dbar plot = special legend
        plt.legend(legendArtistProxyShapes,
                options.custom_legend_labels if options.custom_legend_labels != None else tagLabels)
    elif options.custom_legend_labels != None and len(options.custom_legend_labels) > 1 and not skipLegend:
        if options.horizontal:
            handles, labels = plt.gca().get_legend_handles_labels()
            plt.legend(reversed(handles), reversed(labels))
        else:
            plt.legend(loc = 'upper right')

    if options.horizontal:
        plt.ylabel(options.xlabel)
        plt.xlabel(options.ylabel)
    else:
        plt.xlabel(options.xlabel)
        plt.ylabel(options.ylabel)
    plt.title(options.title)






    myax = plt.gca()
    #mysubax = zoomed_inset_axes(myax, 3, loc=2, bbox_to_anchor=(60,100,10,10), borderpad=2) # zoom = 6
    mysubax = zoomed_inset_axes(myax, 12, bbox_to_anchor=(70,160), loc=2, borderpad=0) # zoom = 6
    if options.title == "Network IO per Iteration":
      mysubax = zoomed_inset_axes(myax, 12, bbox_to_anchor=(64,138), loc=2, borderpad=0) # zoom = 6

    for tagIndex, t in enumerate(tags):
        print("plotting tag: " + t.name)

        if t.tagType == TAG_TYPE.XY:
            exes = []
            whys = []
            for row in data:
                if row[3] == t.name and (options.xfilter == None or row[4] in options.xfilter):
                    exes.append(float(row[4]));
                    whys.append(float(row[5]));
            exes, whys = (np.array(t) for t in zip(*sorted(zip(exes, whys))))
            dedup_exes = np.unique(exes)
            dedup_whys = np.empty(dedup_exes.shape)
            for i, x in enumerate(dedup_exes):
                dedup_whys[i] = np.mean(whys[exes == x]) # / 1e9
            if options.title == "Network IO per Iteration":
              dedup_whys /= 1e9

            if t.name == 'emp_float_lm_time_vs_points_per_loc_itr_latency':
                edgecolor=dracula.COLORS[tagIndex+options.color_offset]
                import colorsys
                h, l, s = colorsys.rgb_to_hls(*matplotlib.colors.to_rgb(edgecolor))
                edgehighlight = colorsys.hls_to_rgb(h, min(1, l * 0.5), s = s)

                plt.plot(dedup_exes, dedup_whys,
                    label=(options.custom_legend_labels[tagIndex] if options.custom_legend_labels != None else t.name),
                    color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None),
                    linestyle = '--', lw=2, path_effects=[pe.Stroke(linewidth=3, foreground=edgehighlight), pe.Normal()])
            else:
                plt.plot(dedup_exes, dedup_whys,
                    color=(dracula.COLORS[tagIndex+options.color_offset] if options.color_theme == "dracula" else None))

    mysubax.set_xlim(6, 12)
    mysubax.set_xticks([6, 9, 12])
    plt.xticks(visible=True)
    plt.yticks(visible=True)
    if options.title == "Network IO per Iteration":
        mysubax.set_ylim(0, .000003)
        #mysubax.set_yticks([1, 2])
    else:
        mysubax.set_ylim(1, 2)
        mysubax.set_yticks([1, 2])
    #mark_inset(myax, mysubax, loc1=2, loc2=3, fc="none", ec="0.5")
    mark_inset(myax, mysubax, loc1=2, loc2=3, fc="none", ec="0.7")

    from matplotlib.ticker import StrMethodFormatter
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places








    plt.tight_layout()
    if options.show:
        plt.show()
    if options.log_scale:
        if options.horizontal:
            plt.gca().set_xscale('log')
        else:
            plt.gca().set_yscale('log')
    plt.savefig(options.graphpath)
    print("saved figure to " + options.graphpath)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="JLog plotter")
    parser.add_argument('--csvlog',
                    required=True,
                    nargs="+",
                    help=("Path to log file."))
    parser.add_argument('--graphpath',
                    required=True,
                    help=("Path to store output graphs."))
    parser.add_argument('--title',
                    required=False,
                    default="",
                    help=("Graph Title."))
    parser.add_argument('--xlabel',
                    required=False,
                    default="",
                    help=("Label for plot x axis."))
    parser.add_argument('--ylabel',
                    required=False,
                    default="",
                    help=("Label for plot y axis."))
    parser.add_argument('--zlabel',
                    required=False,
                    default="",
                    help=("Label for plot z axis."))
    parser.add_argument('--only-tags',
                    required=False,
                    nargs="+",
                    help=("Only plot specific (space seperated) tags in csv log files. Default plots all tags."))
    parser.add_argument('--custom-legend-labels',
                    required=False,
                    nargs="+",
                    help=("Attach a custom legend label to each tag's plot."))
    parser.add_argument('--color-theme',
                    required=False,
                    help=("Use a custom color theme - dracula suported"))
    parser.add_argument('--color-offset',
                    required=False,
                    type=int,
                    default=0,
                    help=("Offset starting color in color theme"))
    parser.add_argument('--projection',
                    required=False,
                    type=int,
                    help=("Index of x/y projection for 3d to 2d projections"))
    parser.add_argument('--azimuth',
                    required=False,
                    type=int,
                    help=("Azimuth sets view for 3D plot"))
    parser.add_argument('--elevation',
                    required=False,
                    type=int,
                    help=("Elevation sets view for 3D plot"))
    parser.add_argument('--fig-w',
                    required=False,
                    type=float,
                    default=5,
                    help=("Override figure width"))
    parser.add_argument('--fig-h',
                    required=False,
                    type=float,
                    default=4,
                    help=("Override figure height"))
    parser.add_argument('--xfilter',
                    required=False,
                    nargs="+",
                    help=("Only use data which match one of these x values."))
    parser.add_argument('--rotate_x_labels',
                    required=False,
                    action="store_true",
                    help=("Rotate x axis labels."))
    parser.add_argument('--hide-legend',
                    required=False,
                    action="store_true",
                    help=("Hide legend."))
    parser.add_argument('--bar-label-dict',
                    required=False,
                    type=json.loads,
                    help=("Add labels to bars e.g. '{\"tagName\": [\"l1\", \"l2\", \"l3\"]}' "))
    parser.add_argument('--hatching',
                    required=False,
                    help=("Add hatching to bars e.g. \"||||\""))
    parser.add_argument('--show',
                    required=False,
                    action="store_true",
                    help=("Show plots as they are generated."))
    parser.add_argument('--print-raw',
                    required=False,
                    action="store_true",
                    help=("Print raw plotted data to command line."))
    parser.add_argument('--log-scale',
                    required=False,
                    action="store_true",
                    help=("Log scale."))
    parser.add_argument('--horizontal',
                    required=False,
                    action="store_true",
                    help=("Plot bar graphs horizontally."))
    options = parser.parse_args()

    try:
        csvFilePaths = options.csvlog
        for p in csvFilePaths:
            with open(p, "r") as f:
                print("found file " + p)
    except FileNotFoundError:
        print("file not found " + p)
        exit()

    if options.only_tags != None:
        print("Only plotting specific tags: ")
        print(options.only_tags)

    tags = getTags(csvFilePaths, options.only_tags)

    if tags == None or len(tags) == 0:
        print("ERROR - No tags found - Nothing to plot")
        exit()

    if options.custom_legend_labels != None:
        if len(options.custom_legend_labels) == len(tags):
            print("Custom Tag Labels Enabled")
        else:
            print("Wrong Number of Custom Tag Labels")
            print(options.custom_legend_labels)
            print(tags)
            exit()

    plot(csvFilePaths, options, tags)

    if os.path.exists(tmpFile):
        os.remove(tmpFile)
