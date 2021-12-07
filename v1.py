#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program will open a file whose name is specified by the user and
convert each column into a list which can then be filtered by
iterating through the list.

Currently supported filters:
    Rolling Average
    Savitzky Golay
    CFC 180*
    CFC 600*
    CFC 1000*
    Custom CFC Filter*
    Cropping of data with no filtering

    *I believe the CFC filters are put through 4th order butterworth filters; see lines 177 - 190
    *signal.butter(2) generates second-order sections representation of IIR filter
    *and signal.sosfiltfilt cascades two 2-pole filters forwards and then backwards which I believe accomplishes
    *4th order requirements

    *Further information on these functions can be found at https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfiltfilt.html
    *and https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html


The first column should be the acceleration data in g's and the second
should be force in Newtons.

"""
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.signal import savgol_filter
from scipy import signal
from datetime import datetime
import numpy as np

filtered_t_f = False


######################################################## Deletes data starting from the left ###########################################
def deleteBeginningData(beginning, column1, column2):
    bigTempList = []
    tempList1 = []
    tempList2 = []
    bigTempList.clear()  # making SURE list is cleared, had issues during testing without this
    tempList1.clear()
    tempList2.clear()
    for j in range(beginning, (len(column1) - 1)):
        tempList1.append(column1[j])
        tempList2.append(column2[j])
    bigTempList.append(tempList1)
    bigTempList.append(tempList2)
    return bigTempList


######################################################## Deletes data starting from the right ##########################################
def deleteEndData(end, column1, column2):
    bigTempList = []
    tempList1 = []
    tempList2 = []
    bigTempList.clear()  # making SURE list is cleared, had issues during testing without this
    tempList1.clear()
    tempList2.clear()
    for n in range(end):
        tempList1.append(column1[n])
        tempList2.append(column2[n])
    bigTempList.append(tempList1)
    bigTempList.append(tempList2)
    return bigTempList


######################################################## Produces plots of unfiltered data ###########################################
def producePlots(plot_x_size, plot_y_size):  # produces plots of unfiltered data during data selection
    plt.figure(figsize=(plot_x_size, plot_y_size))
    plt.subplot(1, 2, 1)
    plt.plot(column1)
    plt.title("Acceleration data")
    plt.xlabel("Time")
    plt.ylabel("Acceleration (g)")
    plt.subplot(1, 2, 2)
    plt.plot(column2)
    plt.title("Force data")
    plt.xlabel("Time")
    plt.ylabel("Force (N)")
    plt.show()


######################################################## Produces plots of filtered data #############################################
def produceFilteredPlots(plot_x_size, plot_y_size):
    plt.figure(figsize=(plot_x_size, plot_y_size))
    plt.subplot(1, 2, 1)
    plt.plot(biglist2[0])
    plt.title("Filtered Acceleration")
    plt.xlabel("Time (microseconds?)")
    plt.ylabel("Acceleration (m/s^2)")
    plt.subplot(1, 2, 2)
    plt.plot(biglist2[1])
    plt.title("Filtered Force")
    plt.xlabel("Time (microseconds?)")
    plt.ylabel("Acceleration (m/s^2)")
    plt.show()


######################################################## Crops unfiltered data by deleting unneeded datapoints ##################################
def dataSelection(col1, col2, filtered_t_f, plot_x_size, plot_y_size):
    global column1
    global column2
    smallerList = []
    cont = input("Would you like to remove more unneeded data? (y/n): ")
    valid = False
    while valid == False:
        if ((cont == 'y') or (cont == 'n')):
            valid = True
        else:
            print("INVALID INPUT\n")
            cont = input("Would you like to remove more unneeded data? (y/n): ")
    while cont == 'y':
        whichSide = input("Enter 'begin' to remove data at the beginning or 'end' to remove data at the end: ")
        smallerList.clear()
        if whichSide == 'begin':  # delete data from left
            beginning = int(input("Enter x-axis value where desired data begins: "))
            smallerList = deleteBeginningData(beginning, column1, column2)  # holds the column's data after cropping
            column1 = smallerList[0]  # column1 and 2 now have cropped data
            column2 = smallerList[1]
            if filtered_t_f == False:
                producePlots(plot_x_size, plot_y_size)
            else:
                produceFilteredPlots(plot_x_size, plot_y_size)

        elif whichSide == 'end':  # delete data from right
            end = int(input("Enter x-axis value where desired data ends: "))
            smallerList = deleteEndData(end, column1, column2)
            column1 = smallerList[0]  # holds the column's data after cropping
            column2 = smallerList[1]  # column1 and 2 now have cropped data
            if filtered_t_f == False:
                producePlots(plot_x_size, plot_y_size)
            else:
                produceFilteredPlots(plot_x_size, plot_y_size)
        else:
            print("INVALID INPUT\n")
        cont = input("Delete more data? (y/n): ")


######################################################## Crops filtered data by deleting unneeded datapoints ##################################
def filteredDataSelection(col1, col2, filtered_t_f, plot_x_size, plot_y_size):
    global biglist2
    cont = 'y'
    cont = input("Would you like to remove more unneeded data? (y/n): ")
    smallerList = []
    while cont == 'y':
        whichSide = input("Enter 'begin' to remove data at the beginning or 'end' to remove data at the end: ")
        smallerList.clear()
        if whichSide == 'begin':  # delete data from left
            beginning = int(input("Enter x-axis value where desired data begins: "))
            smallerList = deleteBeginningData(beginning, biglist2[0],
                                              biglist2[1])  # holds the column's data after cropping
            biglist2[0] = smallerList[0]  # column1 and 2 now have cropped data
            biglist2[1] = smallerList[1]
            if filtered_t_f == False:
                producePlots(plot_x_size, plot_y_size)
            else:
                produceFilteredPlots(plot_x_size, plot_y_size)

        elif whichSide == 'end':  # delete data from right
            end = int(input("Enter x-axis value where desired data ends: "))
            smallerList = deleteEndData(end, biglist2[0], biglist2[1])
            biglist2[0] = smallerList[0]  # holds the column's data after cropping
            biglist2[1] = smallerList[1]  # column1 and 2 now have cropped data
            producePlots(plot_x_size, plot_y_size)
            if filtered_t_f == False:
                producePlots(plot_x_size, plot_y_size)
            else:
                produceFilteredPlots(plot_x_size, plot_y_size)
        else:
            print("INVALID INPUT\n")
        cont = input("Delete more data? (y/n): ")


######################################################## ROLLING AVG FILTER FUNCTION #############################################
def rollingAverage(biglist, windowSize):
    global biglist2
    biglist2 = []
    filteredCol1 = []  # same content as biglist but filtered
    filteredCol2 = []
    numRows = len(biglist[1])  # number of rows after extra data deleted
    iteration = 0
    iterationLimit = (1 + numRows - windowSize)  # max rolling average datapoints that can be generated ...
    # using rolling average filter
    addedWindow1 = 0
    addedWindow2 = 0
    # start = windowSize - 1                                         # index of where first calculation can begin
    while iteration < iterationLimit:
        for x in range(windowSize):
            addedWindow1 += biglist[0][(x + iteration)]  # adding values in window
            addedWindow2 += biglist[1][(x + iteration)]
        averagedWindow1 = (addedWindow1 / windowSize)  # finding average of all values in window
        averagedWindow2 = (addedWindow2 / windowSize)
        filteredCol1.append(averagedWindow1)
        filteredCol2.append(averagedWindow2)
        addedWindow1 = 0  # reset values for new rolling avg calculation
        addedWindow2 = 0
        averagedWindow1 = 0
        averagedWindow2 = 0
        iteration += 1  # increment interation counter
    biglist2.append(filteredCol1)
    biglist2.append(filteredCol2)


######################################################## DEFINING LISTS TO STORE COLUMN DATA ###################################
column1 = []  # acceleration
column2 = []  # force
readFile = input("--Enter name of file to be read: ")
######################################################## GENERATING NEW FILE NAMES ###################################
now = datetime.now()  # retrieves data and time
dt_string = now.strftime('%m-%d__%H-%M')  # adds date and time to new file name
newFileName = readFile.replace('.csv', '_filtered_')
# newFileName = readFile.replace('.tsv', '_filtered_')
plotnamebase = newFileName + dt_string  # will be used to generate png saved file names
writeFile = plotnamebase + '.csv'  # name of file that will be written
print("\nWrite file name will be ", writeFile, '\n')
#####################################################################################################################################
filterType = int(input(
    "--Enter number of desired filter\n\n1. Rolling Average \n2. Savitzky Golay \n3. CFC 180 \n4. CFC 600 \n5. CFC 1000 \n6. CFC filter with custom cutoff frequency\n7. No filter, data cropping only:  \n"))
validEntry = False
validNum = False
######################################################## MAKING SURE FILTER EXISTS #############################################
while validEntry == False:
    if filterType == 1:
        windowSize = int(input("--Enter window size for moving average: "))  # size of moving average window
        validEntry = True

    elif filterType == 2:
        validEntry = True
        while validNum == False:
            windowSize = int(input(
                "--Enter window size for filtering\n--Must be odd number, entering 0 means only smoothing: "))  # size of moving average window
            if (windowSize % 2) == 0:
                print("ERROR")
                print("WINDOW MUST BE ODD.\n")
            order = int(input("--Enter order for filter\n--Must be less than window size: "))
            if (order >= windowSize):
                print("ERROR")
                print("ORDER MUST BE LESS THAN WINDOW SIZE.\n")
            else:
                validNum = True

    elif ((filterType == 3) or (filterType == 4) or (filterType == 5)):
        validEntry = True
        f = int(input("Enter data acquisition rate (f = 100000 in Sample_data.csv) :"))
        # frequency variable here

    elif filterType == 6:
        validEntry = True
        customFreq = int(input("Enter desired cutoff frequency: "))
        f = int(input("Enter data acquisition rate (f = 100000 in Sample_data.csv) :"))

    elif filterType == 7:
        validEntry = True

    else:
        print("\nERROR: INVALID FILTER SELECTED")
        filterType = int(input(
            "--Enter number of desired filter\n\nCurrently Supported: \n1. Rolling Average \n2. Savitzky Golay \n3. CFC 180 \n4. CFC 600 \n5. CFC 1000 \n6. CFC filter with custom cutoff frequency\n7. No filter, data cropping only:  \n"))
        validEntry = False
######################################################## OPENING CSV AND PUTTING COLUMN DATA INTO LISTS ########################
csv_or_tsv = readFile.find(".tsv")  # determining file type, returns -1 if ".tsv" is not typed
with open(readFile, 'r') as readFile:
    print("Reading File...\n")
    for eachline in readFile:  # converting columns to lists
        if (csv_or_tsv == -1):
            parts = eachline.strip('\n').split(',')  # if file is a csv
        else:
            parts = eachline.strip('\n').split('\t')  # if file is tsv
        column1.append(parts[0])
        column2.append(parts[1])
    del column1[0]  # deleting the header of the column
    del column2[0]
    for i in range(len(column1)):  # converting strings to int/float
        column1[i] = float(column1[i]) * 9.81  # from g to m/s^2
        column2[i] = float(column2[i])
    ######################################################## USER DEFINES PLOT SIZES ###############################################
    plot_x_size = 20
    plot_y_size = 6.5
    default = input("Use default plot sizes (20 x 6.5)? (y/n) ")
    if (default == 'n'):
        plot_x_size = float(input("Enter x value for plot size: "))
        plot_y_size = float(input("Enter y value for plot size: "))
    ######################################################## USER DEFINES AREA OF DESIRED DATA HERE #################################
    producePlots(plot_x_size, plot_y_size)
    beginning = int(input("Enter x-axis value where desired data begins: "))
    smallerList = deleteBeginningData(beginning, column1, column2)
    column1 = smallerList[0]
    column2 = smallerList[1]
    producePlots(plot_x_size, plot_y_size)
    end = int(input("Enter x-axis value where desired data ends: "))
    smallerList = deleteEndData(end, column1, column2)
    column1 = smallerList[0]
    column2 = smallerList[1]
    producePlots(plot_x_size, plot_y_size)
    dataSelection(column1, column2, filtered_t_f, plot_x_size, plot_y_size)
    # del column1[0]                                  # the first row is all 0's so delete it
    # del column2[0]
    print("Creating big list...\n")
biglist = []  # holds both columns of data
biglist.append(column1)
biglist.append(column2)
######################################################## RUN ROLLING AVG FILTER ################################################
if filterType == 1:
    rollingAverage(biglist, windowSize)
    filtered_t_f = True

######################################################## RUN SAVITZKY GOLAY FILTER #############################################
elif filterType == 2:
    biglist2 = []
    savGolFilter1 = savgol_filter(biglist[0], windowSize, order)  # filtered acceleration data
    savGolFilter2 = savgol_filter(biglist[1], windowSize, order)  # filtered force data
    savGolFilter1 = savGolFilter1.tolist()
    savGolFilter2 = savGolFilter2.tolist()
    biglist2.append(savGolFilter1)
    biglist2.append(savGolFilter2)
    filtered_t_f = True

######################################################## RUN CFC 180 FILTER ###################################################
elif filterType == 3:
    oneEighty = signal.butter(2, 300, btype='low', analog=False, output='sos',
                              fs=f)  # (order, cutoff frequency, lowpass, digital filtering, output second order sections, sampling frequency)
    buttFilter1 = signal.sosfiltfilt(oneEighty, biglist[
        0])  # filter once forwards and once backwards to create fourth order with no phase shift
    buttFilter2 = signal.sosfiltfilt(oneEighty, biglist[1])
    biglist2 = []
    buttFilter1 = buttFilter1.tolist()
    buttFilter2 = buttFilter2.tolist()
    biglist2.append(buttFilter1)
    biglist2.append(buttFilter2)
    filtered_t_f = True

######################################################## RUN CFC 600 FILTER ##################################################
elif filterType == 4:
    oneEighty = signal.butter(2, 1000, btype='low', analog=False, output='sos', fs=f)
    buttFilter1 = signal.sosfiltfilt(oneEighty, biglist[0])
    buttFilter2 = signal.sosfiltfilt(oneEighty, biglist[1])
    biglist2 = []
    buttFilter1 = buttFilter1.tolist()
    buttFilter2 = buttFilter2.tolist()
    biglist2.append(buttFilter1)
    biglist2.append(buttFilter2)
    filtered_t_f = True

######################################################## RUN CFC 1000 FILTER ################################################
elif filterType == 5:
    oneEighty = signal.butter(2, 1650, btype='low', analog=False, output='sos', fs=f)
    buttFilter1 = signal.sosfiltfilt(oneEighty, biglist[0])
    buttFilter2 = signal.sosfiltfilt(oneEighty, biglist[1])
    biglist2 = []
    buttFilter1 = buttFilter1.tolist()
    buttFilter2 = buttFilter2.tolist()
    biglist2.append(buttFilter1)
    biglist2.append(buttFilter2)
    filtered_t_f = True

######################################################## RUN CUSTOM CFC FILTER ###############################################
elif filterType == 6:
    oneEighty = signal.butter(2, customFreq, btype='low', analog=False, output='sos', fs=f)
    buttFilter1 = signal.sosfiltfilt(oneEighty, biglist[0])
    buttFilter2 = signal.sosfiltfilt(oneEighty, biglist[1])
    # indexCol = [ ]
    biglist2 = []
    # for z in range(len(buttFilter1)):                              # creating index column (left code here in case needed in future)
    #   indexCol.append(z)
    # biglist2.append(indexCol)
    buttFilter1 = buttFilter1.tolist()
    buttFilter2 = buttFilter2.tolist()
    biglist2.append(buttFilter1)
    biglist2.append(buttFilter2)
    filtered_t_f = True

######################################################## NO FILTERING ########################################################
elif filterType == 7:
    biglist2 = []
    biglist2.append(biglist[0])
    biglist2.append(biglist[1])
produceFilteredPlots(plot_x_size, plot_y_size)
filteredDataSelection(biglist2[0], biglist2[1], filtered_t_f, plot_x_size, plot_y_size)
velocity = integrate.cumtrapz(biglist2[0], dx=1.0)  # calculating velocity by taking integral of acceleration
npdisplacement = integrate.cumtrapz(velocity, dx=1.0)  # calculating displacement by taking integral of velocity
biglist2[0] = biglist2[0][0: (len(biglist2[0]) - 2)]  # making lists the same size for csv writing, loses 2 datapoints
biglist2[1] = biglist2[1][0: (len(biglist2[0]))]
velocity = velocity[0: (len(velocity) - 1)]  # making lists the same size for csv writing, loses 1 datapoint
displacement = []
for x in range(len(npdisplacement)):
    num = float(npdisplacement[x])
    displacement.append(num)
combinedList = np.append(displacement, biglist2[1])
midpoint = (len(combinedList) // 2)
force = combinedList[0: midpoint]
disp = combinedList[midpoint:]
function = np.polyfit(force, disp, 2)  # finds 2nd order polynomial coefficients for force v displacement plot
totalEnergyAbsorbed = ((function[0] * (force[-1]) ** 2) + function[1] * force[-1] + function[2]) - (
        (function[0] * (force[0]) ** 2) + function[1] * force[0] + function[2])
# plugging in last x value into polynomial minus first x value plugged into polynomial
print("Total energy absorbed: ", totalEnergyAbsorbed)
finalList = []
finalList.append(biglist2[0])
finalList.append(velocity)
finalList.append(disp)
finalList.append(biglist2[1])
######################################################## WRITING FILTERED DATA TO NEW FILE ###################################
with open(writeFile, 'w') as writeFile:
    print("Writing filtered data file...\n")
    writeFile.write('Acceleration (m/s^2),')
    writeFile.write('Velocity (m/s),')
    writeFile.write('Displacement (m),')
    writeFile.write('Force (N),')
    writeFile.write('\n')
    for x in range(len(finalList[1])):  # will iterate as many times as there are rows in csv
        for y in range(len(finalList)):  # will iterate as many times as there are columns in csv (4)
            if y == 3:
                writeFile.write(str(finalList[y][x]))

            else:
                writeFile.write(str(finalList[y][x]) + ',')
        writeFile.write('\n')
######################################################## FINDING CRITICAL VALUES AND PUTTING INTO CSV ##############################
minAcceleration = min(biglist2[0])
maxAcceleration = max(biglist2[0])
minForce = min(biglist2[1])
maxForce = max(biglist2[1])
minVelocity = min(velocity)
maxVelocity = max(velocity)
minDisplacement = min(displacement)
maxDisplacement = max(displacement)
critvalues = []
firstcol = ['\0', 'Acceleration (m/s^2)', 'Force (N)', 'Velocity (m/s)', 'Displacement (m)']
secondcol = ['Minimum', minAcceleration, minForce, minVelocity, minDisplacement]
thirdcol = ['Maximum', maxAcceleration, maxForce, maxVelocity, maxDisplacement]
fourthcol = ['Total Energy Absorbed', totalEnergyAbsorbed, '\0', '\0', '\0']
critvalues.append(firstcol)
critvalues.append(secondcol)
critvalues.append(thirdcol)
critvalues.append(fourthcol)
critFileName = newFileName + dt_string + '_criticalvalues.csv'
with open(critFileName, 'w') as writeFile:
    print("Writing critical values file...\n")
    for x in range(len(critvalues[1])):  # will iterate as many times as there are rows in csv
        for y in range(len(critvalues)):  # will iterate as many times as there are columns in csv (4)
            if y == 3:
                writeFile.write(str(critvalues[y][x]))

            else:
                writeFile.write(str(critvalues[y][x]) + ',')
        writeFile.write('\n')
print("Program Complete!\n")
######################################################## SAVING AND DISPLAYING PLOTS ########################################################
plt.figure(figsize=(plot_x_size, plot_y_size))
plt.plot(biglist2[0])
plt.title("Filtered Acceleration v Time")
plt.xlabel("Time")
plt.ylabel("Acceleration (g)")
saveName = plotnamebase + '_filteredAcceleration.png'
plt.savefig(saveName)
plt.figure(figsize=(plot_x_size, plot_y_size))
plt.plot(biglist2[1])
plt.title("Filtered Force v Time")
plt.xlabel("Time")
plt.ylabel("Force (N)")
saveName = plotnamebase + '_filteredForce.png'
plt.savefig(saveName)
plt.figure(figsize=(plot_x_size, plot_y_size))
plt.plot(velocity)
plt.title("Velocity v Time")
plt.xlabel("Time")
plt.ylabel("Velocity (N)")
saveName = plotnamebase + '_velocity.png'
plt.savefig(saveName)
plt.figure(figsize=(plot_x_size, plot_y_size))
plt.plot(displacement)
plt.title("Displacement v Time")
plt.xlabel("Time")
plt.ylabel("Displacement (m)")
saveName = plotnamebase + '_displacement.png'
plt.savefig(saveName)
plt.figure(figsize=(plot_x_size, plot_y_size))
plt.plot(force, disp)
plt.title("Force vs Displacemnt")
plt.xlabel("Displacement (m)")
plt.ylabel("Force (N)")
saveName = plotnamebase + '_forceVdisplacement.png'
plt.savefig(saveName)
print("Filtered data written to", plotnamebase)
if filterType == 1:
    print("Filter used was Rolling Average")
elif filterType == 2:
    print("Filter used was Savitzky Golay")
elif filterType == 3:
    print("Filter used was CFC 180")
elif filterType == 4:
    print("Filter used was CFC 600")
elif filterType == 5:
    print("Filter used was CFC 1000")
elif filterType == 6:
    print("Filter used was Custom CFC with freq cutoff ", customFreq, "Hz")
elif filterType == 7:
    print("Data not filtered, only cropped")
