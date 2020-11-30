import glob
import os
from datetime import datetime
import numpy as np
from nptdms import TdmsFile
from scipy.signal import savgol_filter
import pandas as pd
from scipy.stats import stats
import matplotlib.pyplot as plt


def smooth_fun(arr, span=20):
    smooth_batch = int(len(arr) / span)  # find smoothing window (ex. every 5 points)
    if smooth_batch % 2 == 0:  # smoothing window must be odd
        smooth_batch += 1
    return savgol_filter(arr, smooth_batch, 2)


loading_color = '#2039aa'
unloading_color = '#9c2a6f'


class Curve:
    def __init__(self, position=np.array([]), load=np.array([]), position_discharge=np.array([]),
                 load_discharge=np.array([])):
        self.load = load
        self.load_discharge = load_discharge
        self.position = position
        self.position_discharge = position_discharge

    def interpolate(self, loading_len, unloading_len, cutoff=False):
        load = smooth_fun(self.load)
        position = smooth_fun(self.position)
        load_discharge = smooth_fun(self.load_discharge)
        position_discharge = smooth_fun(self.position_discharge)

        if cutoff:
            position = np.flip(position)
            load = np.flip(load)
        else:
            position_discharge = np.flip(position_discharge)
            load_discharge = np.flip(load_discharge)
        x_vals = np.linspace(min(self.position), max(self.position), loading_len)
        x_vals_discharge = np.linspace(min(self.position_discharge), max(self.position_discharge), unloading_len)

        y_interp = np.interp(x_vals, position, load)
        y_interp_discharge = np.interp(x_vals_discharge, position_discharge, load_discharge)

        return Curve(load=y_interp, load_discharge=y_interp_discharge,
                     position=x_vals, position_discharge=x_vals_discharge)

    def plot(self):
        plt.xlabel('Load [grF]')
        plt.ylabel('Displacement [μm]\nDcorrected')
        plt.plot(self.load, self.position, color=loading_color, label="loading")
        plt.plot(self.load_discharge, self.position_discharge,
                 color=unloading_color,
                 label="unloading")

        plt.legend(bbox_to_anchor=(.7, 1), loc='upper left')
        plt.show()


class Experiment:

    def __init__(self, path, derivative_offset=-0.1):
        try:
            self.failed = False
            self.status = "UNKNOWN"
            self.cut_off_position = 0

            self.filename = os.path.basename(path)
            tdms_file = TdmsFile.read(path)
            group = tdms_file['Untitled']
            position = group['Position'][:]
            load = group['Load'][:]
            peek_pos_index = np.argmin(position) + 1  # for discharge

            self.original = Curve(load=load[:peek_pos_index + 1],
                                  load_discharge=load[peek_pos_index:],
                                  position=position[:peek_pos_index + 1],
                                  position_discharge=position[peek_pos_index:])

            self.derivative_offset = derivative_offset
            self.cut_off_curve = self.cut_off_inactive()
            self.mirror_offset_curve = self.mirror_and_offset()
            self.stiffness_curve = self.stiffness()
            self.validity_check()
        except:
            self.failed = True

    def process(self):
        if not self.failed:
            self.cut_off_curve = self.cut_off_inactive()
            self.mirror_offset_curve = self.mirror_and_offset()
            self.stiffness_curve = self.stiffness()

    def cut_off_inactive(self):
        if self.cut_off_position != 0:  # Manual cutoff
            tmp_pos = self.original.position[self.original.position > self.cut_off_position]

            if len(tmp_pos) > 0:  # if not empty
                offset_position = np.min(tmp_pos)  # value of the discharge position

                offset_index_position = \
                    np.where(self.original.position == offset_position)[0][0]
            else:
                offset_index_position = 0

            tmp_pos = self.original.position_discharge[self.original.position_discharge > self.cut_off_position]
            if len(tmp_pos) > 0:  # if not empty
                offset_position_discharge = np.min(tmp_pos)  # value of the discharge position
                offset_index_position_discharge = \
                    np.where(self.original.position_discharge == offset_position_discharge)[0][0]
            else:
                offset_index_position_discharge = -1

            return Curve(load=self.original.load[offset_index_position:],
                         load_discharge=self.original.load_discharge[:offset_index_position_discharge],
                         position=self.original.position[offset_index_position:],
                         position_discharge=self.original.position_discharge[:offset_index_position_discharge])

        smooth_load = smooth_fun(self.original.load)
        smooth_position = smooth_fun(self.original.position)

        dl_dp = np.gradient(smooth_load, smooth_position)  # get the derivative

        smooth_dl_dp = smooth_fun(dl_dp)
        # find the closest value to derivative_offset
        offset_value = np.min(smooth_dl_dp[smooth_dl_dp > self.derivative_offset])
        # find index of this value, aka the load offset
        offset_index = np.where(smooth_dl_dp == offset_value)[0][0]

        offset_position = self.original.position[offset_index]  # value of the charge position
        offset_position_discharge = np.min(
            self.original.position_discharge[
                self.original.position_discharge > offset_position])  # value of the discharge position
        offset_index_position_discharge = np.where(self.original.position_discharge == offset_position_discharge)[0][0]

        # keep only the values after the offset, both on load and position
        # for the discharge we offset based on the position we calculate
        return Curve(load=self.original.load[offset_index:],
                     load_discharge=self.original.load_discharge[:offset_index_position_discharge],
                     position=self.original.position[offset_index:],
                     position_discharge=self.original.position_discharge[:offset_index_position_discharge])

    def mirror_and_offset(self):
        n0 = self.cut_off_curve.position[0]

        return Curve(load=self.cut_off_curve.load,
                     load_discharge=self.cut_off_curve.load_discharge,
                     position=n0 - self.cut_off_curve.position,
                     position_discharge=n0 - self.cut_off_curve.position_discharge)

    def stiffness(self):
        return Curve(load=self.mirror_offset_curve.load,
                     load_discharge=self.mirror_offset_curve.load_discharge,
                     position=self.mirror_offset_curve.position -
                              self.mirror_offset_curve.load * 1.757,
                     position_discharge=self.mirror_offset_curve.position_discharge -
                                        self.mirror_offset_curve.load_discharge * 1.757)

    def validity_check(self):
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.stiffness_curve.load,
                                                                       self.stiffness_curve.position)
        self.status = "TRASH" if slope < 0 else "KEEP"

    def plot_cutoff(self, axes, loading=True, unloading=True):
        axes.set_title('Initial Measurement')
        axes.set_xlabel('Displacement [μm]')
        axes.set_ylabel('Load [grF]')
        if loading:
            axes.plot(self.original.position, self.original.load, color='#6a84fb', label="loading offset")
            axes.plot(self.cut_off_curve.position, self.cut_off_curve.load, color=loading_color, label="loading")

        if unloading:
            axes.plot(self.original.position_discharge, self.original.load_discharge, color='#e460b0',
                      label="uloading offset")
            axes.plot(self.cut_off_curve.position_discharge, self.cut_off_curve.load_discharge, color=unloading_color,
                      label="unloading")

        axes.text(0, 1.03, self.status, transform=axes.transAxes,
                  bbox=dict(facecolor=self.status_color(), alpha=.5))

        axes.legend(bbox_to_anchor=(.7, 1), loc='upper left', borderaxespad=0.)

    def plot_stiffness(self, axes, loading=True, unloading=True):
        axes.set_title('Stiffness Corrected Measurement')
        axes.set_xlabel('Load [grF]')
        axes.set_ylabel('Displacement [μm]\nDcorrected')
        if loading:
            axes.plot(self.stiffness_curve.load, self.stiffness_curve.position, color=loading_color, label="loading")
        if unloading:
            axes.plot(self.stiffness_curve.load_discharge, self.stiffness_curve.position_discharge,
                      color=unloading_color,
                      label="unloading")
        axes.text(0, 1.03, self.status, transform=axes.transAxes,
                  bbox=dict(facecolor=self.status_color(), alpha=.5))
        axes.legend(bbox_to_anchor=(.7, 1), loc='upper left', borderaxespad=0.)

    def status_color(self):
        if self.status == "UNKNOWN":
            return '#a9a9a9'
        elif self.status == "KEEP":
            return '#69ffa1'
        else:
            return '#ff6969'


class Experiments:
    def __init__(self, dataset_folder):
        self.dataset_folder = dataset_folder
        self.datasets = list(glob.glob(self.dataset_folder + "/*.tdms"))
        self.experiments = list(map(lambda path: Experiment(path=path), self.datasets))

        self.pos = 0

        self.mean_cut_off_curve = Curve()
        self.mean_stiffness_curve = Curve()

        # self.mean_cut_off()
        self.find_mean_of_all()

    def mean(self, the_curve, experiments):
        arrs = list(map(the_curve, experiments))
        return self.avgNestedLists(arrs)

    def avgNestedLists(self, nested_vals):
        """
        Averages a 2-D array and returns a 1-D array of all of the columns
        averaged together, regardless of their dimensions.
        """
        output = []
        maximum = 0
        for lst in nested_vals:
            if len(lst) > maximum:
                maximum = len(lst)
        for index in range(maximum):  # Go through each index of longest list
            temp = []
            for lst in nested_vals:  # Go through each list
                if index < len(lst):  # If not an index error
                    temp.append(lst[index])
            output.append(np.nanmean(temp))
        return output

    def find_max_len(self, curve, experiments):
        data = list(map(curve, experiments))
        return max(list(map(lambda x: len(x), data)))

    def movingaverage(self, interval, window_size):
        window = np.ones(int(window_size)) / float(window_size)
        return np.convolve(interval, window, 'same')

    def find_mean_of_all(self):
        experiments = list(filter(lambda x: not x.failed and x.status != "TRASH", self.experiments))
        loading_len = self.find_max_len(lambda x: x.cut_off_curve.position, experiments)
        unloading_len = self.find_max_len(lambda x: x.cut_off_curve.position_discharge, experiments)
        experiments_cut_off_interp = list(
            map(lambda x: x.cut_off_curve.interpolate(loading_len, unloading_len, cutoff=True), experiments))
        experiments_stiffness_interp = list(
            map(lambda x: x.stiffness_curve.interpolate(loading_len, unloading_len, cutoff=False), experiments))
        print(experiments_stiffness_interp)

        cut_off_loads = list(map(lambda x: x.load, experiments_cut_off_interp))
        cut_off_loads_discharge = list(map(lambda x: x.load_discharge, experiments_cut_off_interp))
        cut_off_positions = list(map(lambda x: x.position, experiments_cut_off_interp))
        cut_off_positions_discharge = list(map(lambda x: x.position_discharge, experiments_cut_off_interp))

        stiffness_loads = list(map(lambda x: x.load, experiments_stiffness_interp))
        stiffness_loads_discharge = list(map(lambda x: x.load_discharge, experiments_stiffness_interp))
        stiffness_positions = list(map(lambda x: x.position, experiments_stiffness_interp))
        stiffness_positions_discharge = list(map(lambda x: x.position_discharge, experiments_stiffness_interp))

        self.mean_cut_off_curve = Curve(
            load=np.average(cut_off_loads, axis=0),
            load_discharge=np.average(cut_off_loads_discharge, axis=0),
            position=np.average(cut_off_positions, axis=0),
            position_discharge=np.average(cut_off_positions_discharge, axis=0))

        self.mean_stiffness_curve = Curve(
            load=np.average(stiffness_loads, axis=0),
            load_discharge=np.average(stiffness_loads_discharge, axis=0),
            position=np.average(stiffness_positions, axis=0),
            position_discharge=np.average(stiffness_positions_discharge, axis=0))

    def next(self):
        self.pos = min(self.pos + 1, len(self.experiments) - 1)

    def prev(self):
        self.pos = max(self.pos - 1, 0)

    def set_status(self, status):
        self.get_selected().status = status

    def process(self):
        self.get_selected().process()

    def get_selected_min_max(self):
        return np.max(self.get_selected().original.position), \
               np.min(self.get_selected().original.position)

    def get_selected(self):
        return self.experiments[self.pos]

    def set_cut_off_position(self, value):
        self.get_selected().cut_off_position = value

    def get_stats(self):
        n_of_experiments = len(self.experiments)
        accepted = 0
        trashed = 0
        remaining = n_of_experiments
        failed = 0

        for experiment in self.experiments:
            if experiment.failed:
                failed += 1
                remaining -= 1
            else:
                if experiment.status == "TRASH":
                    trashed += 1
                    remaining -= 1
                elif experiment.status == "KEEP":
                    accepted += 1
                    remaining -= 1

        return f'{self.pos + 1}/{n_of_experiments}   |   accepted: {accepted}   |' \
               f'   trashed: {trashed}   |   remaining: {remaining}   |   ' \
               f'failed: {failed}   |   {self.get_selected().filename}'

    def plot_all_stiffness_points(self, axes, loading=True, unloading=True, mean=False, trashed=False):
        for experiment in self.experiments:
            if experiment.failed or (not trashed and experiment.status == "TRASH"):
                continue
            if loading:
                color = "#ff2d2d" if experiment.status == "TRASH" else loading_color
                axes.scatter(experiment.stiffness_curve.load, experiment.stiffness_curve.position, c=color, s=1,
                             alpha=0.3, edgecolors='none')

            if unloading:
                color = "#ffd22d" if experiment.status == "TRASH" else unloading_color
                axes.scatter(experiment.stiffness_curve.load_discharge, experiment.stiffness_curve.position_discharge,
                             c=color, s=1,
                             alpha=0.3, edgecolors='none')
        if mean:
            if loading:
                axes.plot(self.mean_stiffness_curve.load, self.mean_stiffness_curve.position,
                          color='green',
                          label="mean loading")
            if unloading:
                axes.plot(self.mean_stiffness_curve.load_discharge, self.mean_stiffness_curve.position_discharge,
                          color='yellow',
                          label="mean unloading")

    def plot_all_cutoff_points(self, axes, loading=True, unloading=True, trashed=False):
        for experiment in self.experiments:
            if experiment.failed or (not trashed and experiment.status == "TRASH"):
                continue
            if loading:
                color = "#ff2d2d" if experiment.status == "TRASH" else loading_color
                axes.scatter(experiment.cut_off_curve.position, experiment.cut_off_curve.load, c=color, s=1,
                             alpha=0.3, edgecolors='none')
            if unloading:
                color = "#ffd22d" if experiment.status == "TRASH" else unloading_color

                axes.scatter(experiment.cut_off_curve.position_discharge, experiment.cut_off_curve.load_discharge,
                             c=color, s=1,
                             alpha=0.3, edgecolors='none')

    def plot_mean_cutoff_points(self, axes, loading=True, unloading=True):
        if loading:
            axes.plot(self.mean_cut_off_curve.position, self.mean_cut_off_curve.load,
                      color='green',
                      label="mean loading")
        if unloading:
            axes.plot(self.mean_cut_off_curve.position_discharge, self.mean_cut_off_curve.load_discharge,
                      color='black',
                      label="mean unloading")

    def plot_mean_stiffness(self, axes, loading=True, unloading=True):
        if unloading:
            axes.plot(self.mean_stiffness_curve.load_discharge, self.mean_stiffness_curve.position_discharge,
                      color='yellow',
                      label="mean unloading")
        if loading:
            axes.plot(self.mean_stiffness_curve.load, self.mean_stiffness_curve.position,
                      color='green',
                      label="mean loading")

    def save_experiments(self):
        timestamp = datetime.now().strftime("%d_%m_%Y__%H_%M_%S")
        path = f"{self.dataset_folder}/results_{timestamp}"
        os.mkdir(path)
        counter = 0
        for experiment in self.experiments:
            if experiment.failed or experiment.status != "KEEP":
                continue
            filename = os.path.splitext(os.path.basename(experiment.filename))[0]
            stiffness_curve = {
                'Load': experiment.stiffness_curve.load,
                'Displacement': experiment.stiffness_curve.position
            }
            stiffness_curve = pd.DataFrame(stiffness_curve, columns=['Load', 'Displacement'])
            stiffness_curve.to_csv(f'{path}/{filename}_loading.csv', index=False, header=True)

            stiffness_curve_discharge = {
                'Load': experiment.stiffness_curve.load_discharge,
                'Displacement': experiment.stiffness_curve.position_discharge
            }

            stiffness_curve = pd.DataFrame(stiffness_curve_discharge, columns=['Load', 'Displacement'])
            stiffness_curve.to_csv(f'{path}/{filename}_unloading.csv', index=False, header=True)
            counter += 1
        return counter, path
