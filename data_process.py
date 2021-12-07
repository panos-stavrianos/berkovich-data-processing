import glob
import os
import traceback
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from nptdms import TdmsFile
from scipy.signal import savgol_filter
from scipy.stats import stats


def sort(a, b):
    c = np.rec.fromarrays([a, b])
    c.sort()
    return c.f0, c.f1


def exp_fit2(x, a, b, c):
    return a * np.exp(-b * x) + c


def exp_fit(x, a, b, c, d, e, f):
    return (a * x) + (b * x ** 2) + (c * x ** 3) + (d * x ** 4) + (e * x ** 5) + f




def smooth_convolve(arr, span=20):
    smooth_batch = int(len(arr) / span)  # find smoothing window (ex. every 5 points)
    if smooth_batch % 2 == 0:  # smoothing window must be odd
        smooth_batch += 1
    if smooth_batch < 2:
        smooth_batch = 3
    return np.convolve(arr, np.ones(smooth_batch) / smooth_batch, mode='valid')



def smooth_fun(arr, span=20):
    smooth_batch = int(len(arr) / span)  # find smoothing window (ex. every 5 points)
    if smooth_batch % 2 == 0:  # smoothing window must be odd
        smooth_batch += 1
    if smooth_batch < 2:
        smooth_batch = 3
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

    def interpolate(self, loading_len, unloading_len):
        load = smooth_fun(self.load)
        position = smooth_fun(self.position)
        load_discharge = smooth_fun(self.load_discharge)
        position_discharge = smooth_fun(self.position_discharge)

        position_discharge = np.flip(position_discharge)
        load_discharge = np.flip(load_discharge)

        load_vals = np.linspace(min(self.load), max(self.load), loading_len)
        load_vals_discharge = np.linspace(min(self.load_discharge), max(self.load_discharge), unloading_len)

        pos_interp = np.interp(load_vals, load, position)
        pos_interp_discharge = np.interp(load_vals_discharge, load_discharge, position_discharge)
        return Curve(load=load_vals, load_discharge=load_vals_discharge,
                     position=pos_interp, position_discharge=pos_interp_discharge)

    def plot(self):
        plt.xlabel('Load [grF]')
        plt.ylabel('Displacement [μm]\nDcorrected')
        plt.xlim([-1, 1.8])
        plt.ylim([20, 30])

        plt.plot(self.load, self.position, color=loading_color, label="loading")
        plt.plot(self.load_discharge, self.position_discharge,
                 color=unloading_color,
                 label="unloading")

        plt.legend(bbox_to_anchor=(.7, 1), loc='upper left')
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()

    def smooth(self, window):
        return Curve(load=smooth_fun(self.load, window), load_discharge=smooth_fun(self.load_discharge, window),
                     position=smooth_fun(self.position, window),
                     position_discharge=smooth_fun(self.position_discharge, window))

    def stiffness(self):
        return Curve(load=self.load,
                     load_discharge=self.load_discharge,
                     position=self.position -
                              self.load * 0.2125,
                     position_discharge=self.position_discharge -
                                        self.load_discharge * 0.2125)

    def offset_load(self, real_max_load=1.5):
        max_load = self.load[-1]
        diff = real_max_load - max_load
        diff = 0
        return Curve(load=self.load + diff,
                     load_discharge=self.load_discharge + diff,
                     position=self.position,
                     position_discharge=self.position_discharge)

    def mirror_and_offset(self):
        n0 = self.position[0]
        return Curve(load=self.load,
                     load_discharge=self.load_discharge,
                     position=n0 - self.position,
                     position_discharge=n0 - self.position_discharge)

    def cut_off_inactive(self, derivative_offset, cut_off_position):
        if cut_off_position != 0:  # Manual cutoff
            tmp_pos = self.position[self.position > cut_off_position]

            if len(tmp_pos) > 0:  # if not empty
                offset_position = np.min(tmp_pos)  # value of the discharge position

                offset_index_position = \
                    np.where(self.position == offset_position)[0][0]
            else:
                offset_index_position = 0

            tmp_pos = self.position_discharge[
                self.position_discharge > cut_off_position]
            if len(tmp_pos) > 0:  # if not empty
                offset_position_discharge = np.min(tmp_pos)  # value of the discharge position
                offset_index_position_discharge = \
                    np.where(self.position_discharge == offset_position_discharge)[0][0]
            else:
                offset_index_position_discharge = -1

            return Curve(load=self.load[offset_index_position:],
                         load_discharge=self.load_discharge[:offset_index_position_discharge],
                         position=self.position[offset_index_position:],
                         position_discharge=self.position_discharge[:offset_index_position_discharge])
        # Automatic cutoff
        dl_dp = np.gradient(self.load, self.position)  # get the derivative

        smooth_dl_dp = smooth_fun(dl_dp)
        # find the closest value to derivative_offset
        offset_value = np.min(smooth_dl_dp[smooth_dl_dp > derivative_offset])
        # find index of this value, aka the load offset
        offset_index = np.where(smooth_dl_dp == offset_value)[0][0]

        offset_position = self.position[offset_index]  # value of the charge position
        offset_position_discharge = np.min(
            self.position_discharge[
                self.position_discharge > offset_position])  # value of the discharge position
        offset_index_position_discharge = \
            np.where(self.position_discharge == offset_position_discharge)[0][0]

        # keep only the values after the offset, both on load and position
        # for the discharge we offset based on the position we calculate
        return Curve(load=self.load[offset_index:],
                     load_discharge=self.load_discharge[:offset_index_position_discharge],
                     position=self.position[offset_index:],
                     position_discharge=self.position_discharge[:offset_index_position_discharge])

    def recalculate_load_upload(self):
        l = np.concatenate([self.load, self.load_discharge])
        p = np.concatenate([self.position, self.position_discharge])
        split_index = np.argmin(p) + 1
        self.load = l[:split_index]
        self.load_discharge = l[split_index:]
        self.position = p[:split_index]
        self.position_discharge = p[split_index:]

    def connect_load_unload(self):
        connected = Curve(load=self.load,
                          load_discharge=self.load_discharge,
                          position=self.position,
                          position_discharge=self.position_discharge)

        connected.load_discharge[0] = self.load[-1]
        connected.position_discharge[0] = self.position[-1]
        return connected


class Experiment:

    def __init__(self, path, derivative_offset=-0.1):
        self.derivative_offset = derivative_offset

        try:
            self.failed = False
            self.status = "UNKNOWN"
            self.cut_off_position = 0

            self.filename = os.path.basename(path)
            tdms_file = TdmsFile.read(path)

            group = tdms_file['Untitled']

            if 'Load' in group:
                load = group['Load'][:]
            else:
                load = group['Untitled 1'][:]

            if 'Position' in group:
                position = group['Position'][:]
            else:
                position = group['Untitled'][:]
            load, position = smooth_convolve(load, 700), smooth_convolve(position, 700)

            peek_load_index = np.argmax(load)  # for discharge
            self.original = Curve(load=load[:peek_load_index],
                                  position=position[:peek_load_index],
                                  load_discharge=load[peek_load_index:],
                                  position_discharge=position[peek_load_index:])
            #self.original.plot()

            load, position = smooth_convolve(load, 70), smooth_convolve(position, 70)

            load, position = smooth_fun(load, 70), smooth_fun(position, 70)
            load, position = smooth_convolve(load, 50), smooth_convolve(position, 50)
            # load, position = smooth_convolve(load, 20), smooth_convolve(position, 20)


            peek_load_index = np.argmax(load)  # for discharge

            self.original_smooth = Curve(load=load[:peek_load_index],
                                  position=position[:peek_load_index],
                                  load_discharge=load[peek_load_index:],
                                  position_discharge=position[peek_load_index:])
            #self.original.plot()
            self.stiffness_curve = None
            self.offset_load = None
            self.cut_off_curve = None
            self.mirror_offset_curve = None

            self.process()
            self.validity_check()
        except Exception:
            print(self.filename)
            traceback.print_exc()
            print('-----------------------------------')
            self.failed = True

    def __repr__(self):
        return f'status: {self.status}'

    def process(self):
        if not self.failed:
            # self.original.plot()

            self.original_stiffness = self.original.stiffness()
            self.stiffness_curve = self.original_smooth.stiffness()

            self.offset_load = self.stiffness_curve.offset_load()
            self.cut_off_curve = self.stiffness_curve.cut_off_inactive(self.derivative_offset, self.cut_off_position)

            self.mirror_offset_curve = self.cut_off_curve.mirror_and_offset()

    def validity_check(self):
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.mirror_offset_curve.load,
                                                                       self.mirror_offset_curve.position)
        self.status = "TRASH" if slope < 0 else "KEEP"

    def plot_cutoff(self, axes, params):
        axes.set_title('Initial Measurement')
        axes.set_xlabel('Displacement [μm]')
        axes.set_ylabel('Load [grF]')
        axes.plot(self.original_stiffness.position, self.original_stiffness.load, color='gray',
                  label="original")
        if params['loading']:
            axes.plot(self.offset_load.position, self.offset_load.load, color='#6a84fb', label="loading offset")
            axes.plot(self.cut_off_curve.position, self.cut_off_curve.load, color=loading_color, label="loading")

        if params['unloading']:
            axes.plot(self.offset_load.position_discharge, self.offset_load.load_discharge, color='#e460b0',
                      label="unloading offset")
            axes.plot(self.cut_off_curve.position_discharge, self.cut_off_curve.load_discharge, color=unloading_color,
                      label="unloading")

        axes.text(0, 1.03, self.status, transform=axes.transAxes,
                  bbox=dict(facecolor=self.status_color(), alpha=.5))

        axes.legend(bbox_to_anchor=(.7, 1), loc='upper left', borderaxespad=0.)

    def plot_stiffness(self, axes, params):
        axes.set_title('Stiffness Corrected Measurement')
        axes.set_xlabel('Load [grF]')
        axes.set_ylabel('Displacement [μm]\nDcorrected')
        if params['loading']:
            axes.plot(self.mirror_offset_curve.load, self.mirror_offset_curve.position, color=loading_color,
                      label="loading")
        if params['unloading']:
            axes.plot(self.mirror_offset_curve.load_discharge, self.mirror_offset_curve.position_discharge,
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
        self.mean_window = 20
        print('find_mean_of_all')
        self.find_mean_of_all()

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

    def find_mean_of_all(self):
        experiments = list(filter(lambda x: not x.failed and x.status != "TRASH", self.experiments))

        loading_len = self.find_max_len(lambda x: x.cut_off_curve.position, experiments)
        unloading_len = self.find_max_len(lambda x: x.cut_off_curve.position_discharge, experiments)
        experiments_stiffness_interp = np.array(list(
            map(lambda x: x.stiffness_curve.interpolate(loading_len, unloading_len), experiments)))

        stiffness_loads = np.array(list(map(lambda x: x.load, experiments_stiffness_interp)))
        stiffness_loads_discharge = list(map(lambda x: x.load_discharge, experiments_stiffness_interp))
        stiffness_positions = np.array(list(map(lambda x: x.position, experiments_stiffness_interp)))
        stiffness_positions_discharge = list(map(lambda x: x.position_discharge, experiments_stiffness_interp))
        self.mean_stiffness_curve = Curve(
            load=np.mean(stiffness_loads, axis=0),
            position=np.mean(stiffness_positions, axis=0),
            load_discharge=np.mean(stiffness_loads_discharge, axis=0),
            position_discharge=np.mean(stiffness_positions_discharge, axis=0))  # .smooth(10)

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

    def get_cut_off_position(self):
        return self.get_selected().cut_off_position

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

    def plot_all_stiffness_points(self, axes, params):
        for experiment in self.experiments:
            if experiment.failed or (not params['trashed'] and experiment.status == "TRASH"):
                continue
            if params['loading']:
                color = "#ff2d2d" if experiment.status == "TRASH" else loading_color
                axes.scatter(experiment.stiffness_curve.load, experiment.stiffness_curve.position, c=color, s=1,
                             alpha=0.3, edgecolors='none')

            if params['unloading']:
                color = "#ffd22d" if experiment.status == "TRASH" else unloading_color
                axes.scatter(experiment.stiffness_curve.load_discharge, experiment.stiffness_curve.position_discharge,
                             c=color, s=1,
                             alpha=0.3, edgecolors='none')

    def plot_all_cutoff_points(self, axes, params):
        for experiment in self.experiments:
            if experiment.failed or (not params['trashed'] and experiment.status == "TRASH"):
                continue
            if params['loading']:
                color = "#ff2d2d" if experiment.status == "TRASH" else loading_color
                axes.scatter(experiment.cut_off_curve.position, experiment.cut_off_curve.load, c=color, s=1,
                             alpha=0.3, edgecolors='none')
            if params['unloading']:
                color = "#ffd22d" if experiment.status == "TRASH" else unloading_color

                axes.scatter(experiment.cut_off_curve.position_discharge, experiment.cut_off_curve.load_discharge,
                             c=color, s=1,
                             alpha=0.3, edgecolors='none')

    def plot_mean_stiffness(self, axes, params):
        self.find_mean_of_all()

        if params['unloading_mean']:
            axes.plot(self.mean_stiffness_curve.load_discharge, self.mean_stiffness_curve.position_discharge,
                      color='yellow',
                      label="mean unloading")
        if params['loading_mean']:
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
            stiffness_curve.to_csv(f'{path}/{filename}_loading.csv', sep=';', index=False, header=True)

            stiffness_curve_discharge = {
                'Load': experiment.stiffness_curve.load_discharge,
                'Displacement': experiment.stiffness_curve.position_discharge
            }

            stiffness_curve = pd.DataFrame(stiffness_curve_discharge, columns=['Load', 'Displacement'])
            stiffness_curve.to_csv(f'{path}/{filename}_unloading.csv', sep=';', index=False, header=True)
            counter += 1

        mean_curve = {
            'Load': self.mean_stiffness_curve.load,
            'Displacement': self.mean_stiffness_curve.position
        }
        mean_curve = pd.DataFrame(mean_curve, columns=['Load', 'Displacement'])
        mean_curve.to_csv(f'{path}/mean_loading.csv', index=False, sep=';', header=True)

        mean_curve_discharge = {
            'Load': self.mean_stiffness_curve.load_discharge,
            'Displacement': self.mean_stiffness_curve.position_discharge
        }
        mean_curve_discharge = pd.DataFrame(mean_curve_discharge, columns=['Load', 'Displacement'])
        mean_curve_discharge.to_csv(f'{path}/mean_unloading.csv', sep=';', index=False, header=True)

        return counter, path
