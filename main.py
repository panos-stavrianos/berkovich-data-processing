import sys
import matplotlib
from PyQt5 import QtWidgets
from PyQt5.QtGui import QKeySequence
from PyQt5.QtWidgets import QFileDialog, QInputDialog, QMessageBox, QAction
from PyQt5.QtWidgets import QLabel

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from data_process import Experiments

matplotlib.use('Qt5Agg')


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=9, height=9, dpi=100):
        experiment = parent.experiments.get_selected()
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        if experiment.failed:
            self.axes.text(0.5, 0.5, f"FAILED\n{experiment.filename}", bbox=dict(facecolor='red', alpha=0.9))
        else:
            if parent.stiffness:
                if parent.all_points:
                    parent.experiments.plot_all_stiffness_points(self.axes, loading=parent.loading,
                                                                 unloading=parent.unloading,
                                                                 trashed=parent.trashed)

                parent.experiments.plot_mean_stiffness(self.axes, loading=parent.loading_mean,
                                                       unloading=parent.unloading_mean)
                experiment.plot_stiffness(self.axes, loading=parent.loading, unloading=parent.unloading)
            else:
                if parent.all_points:
                    parent.experiments.plot_all_cutoff_points(self.axes, loading=parent.loading,
                                                              unloading=parent.unloading,
                                                              trashed=parent.trashed)

                parent.experiments.plot_mean_cutoff_points(self.axes, loading=parent.loading_mean,
                                                           unloading=parent.unloading_mean)
                experiment.plot_cutoff(self.axes, loading=parent.loading, unloading=parent.unloading)

        super(MplCanvas, self).__init__(fig)


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        dialog = QFileDialog()
        self.experiments = Experiments(dialog.getExistingDirectory(self, 'Select the dataset folder'))

        self.stiffness = False
        self.all_points = False
        self.loading = True
        self.unloading = True

        self.trashed = False

        self.loading_mean = False
        self.unloading_mean = False

        self.layout = QtWidgets.QVBoxLayout()
        self.refresh_plot()

        bar = self.menuBar()
        file_menu = bar.addMenu('File')
        actions_menu = bar.addMenu('Actions')
        view_menu = bar.addMenu('View')

        # adding actions to edit menu
        self.create_action('Next', QKeySequence.MoveToNextChar, actions_menu, self.to_next)
        self.create_action('Previous', QKeySequence.MoveToPreviousChar, actions_menu, self.to_prev)
        self.create_action('Set manual offset', QKeySequence("o"), actions_menu, self.set_manual_cutoff)
        self.create_action('Keep', QKeySequence("k"), actions_menu, self.keep_fun)
        self.create_action('Trash', QKeySequence("t"), actions_menu, self.trash_fun)

        self.create_action('Stiffness Corrected Measurement', QKeySequence("TAB"), view_menu, self.toggle_plot_fun,
                           checkable=True)
        self.create_action('All points', QKeySequence("a"), view_menu, self.toggle_all_points_fun,
                           checkable=True)
        self.create_action('Trashed', QKeySequence("d"), view_menu, self.toggle_trashed_fun,
                           checkable=True)
        self.create_action('Loading', QKeySequence("l"), view_menu, self.toggle_loading_fun, checkable=True,
                           checked=True)
        self.create_action('Unloading', QKeySequence("u"), view_menu, self.toggle_unloading_fun,
                           checkable=True, checked=True)
        self.create_action('Loading mean', QKeySequence("Ctrl+l"), view_menu, self.toggle_loading_mean_fun,
                           checkable=True)
        self.create_action('Unloading mean', QKeySequence("Ctrl+u"), view_menu, self.toggle_unloading_mean_fun,
                           checkable=True)

        self.create_action('Open Folder', QKeySequence("Ctrl+O"), file_menu, self.open_folder_fun)
        self.create_action('Save', QKeySequence("Ctrl+S"), file_menu, self.save_fun)
        self.showMaximized()

        self.show()

    def create_action(self, name, shortcut, menu, fun, checkable=False, checked=False):
        action = QAction(name, self, checkable=checkable)
        action.setChecked(checked)

        action.triggered.connect(fun)
        action.setShortcut(shortcut)
        menu.addAction(action)

    def add_status(self):
        label = QLabel(self.experiments.get_stats())
        label.setGeometry(200, 150, 100, 50)
        self.layout.addWidget(label)

        widget = QtWidgets.QWidget()
        widget.setLayout(self.layout)
        self.setCentralWidget(widget)

    def set_manual_cutoff(self):
        max_pos, min_pos = self.experiments.get_selected_min_max()
        d, okPressed = QInputDialog.getDouble(self, "Set Manual Cutoff", "Value: 0 for auto",
                                              max_pos, min_pos, max_pos, 5)
        if okPressed:
            self.experiments.set_cut_off_position(d)
        self.refresh_plot()

    def refresh_plot(self):
        self.clear_layout(self.layout)
        self.experiments.process()
        sc = MplCanvas(self)

        toolbar = NavigationToolbar(sc, self)

        self.layout.addWidget(toolbar)
        self.layout.addWidget(sc)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QtWidgets.QWidget()
        widget.setLayout(self.layout)
        self.setCentralWidget(widget)
        self.add_status()

    def to_next(self):
        print("to_next")
        self.experiments.next()
        self.refresh_plot()

    def to_prev(self):
        print("to_prev")
        self.experiments.prev()
        self.refresh_plot()

    def keep_fun(self):
        print("keep")
        self.experiments.set_status('KEEP')
        self.refresh_plot()

    def trash_fun(self):
        print("trash")
        self.experiments.set_status('TRASH')
        self.refresh_plot()

    def toggle_plot_fun(self):
        print("toggle_plot_fun")
        self.stiffness = not self.stiffness
        self.refresh_plot()

    def toggle_all_points_fun(self):
        print("toggle_plot_fun")
        self.all_points = not self.all_points
        self.refresh_plot()

    def toggle_trashed_fun(self):
        print("toggle_trashed_fun")
        self.trashed = not self.trashed
        self.refresh_plot()

    def toggle_loading_fun(self):
        print("toggle_loading_fun")
        self.loading = not self.loading
        self.refresh_plot()

    def toggle_unloading_fun(self):
        print("toggle_unloading_fun")
        self.unloading = not self.unloading
        self.refresh_plot()

    def toggle_loading_mean_fun(self):
        print("toggle_loading_mean_fun")
        self.loading_mean = not self.loading_mean
        self.refresh_plot()

    def toggle_unloading_mean_fun(self):
        print("toggle_unloading_mean_fun")
        self.unloading_mean = not self.unloading_mean
        self.refresh_plot()

    def toggle_mean_fun(self):
        print("toggle_mean_fun")
        self.mean = not self.mean
        self.refresh_plot()

    def open_folder_fun(self):
        print("open_folder_fun")
        dialog = QFileDialog()
        self.experiments = Experiments(dialog.getExistingDirectory(self, 'Select the dataset folder'))

    def save_fun(self):
        print("save_fun")
        counter, path = self.experiments.save_experiments()
        QMessageBox.about(self, "Saved", f"Αποθηκεύτηκαν {counter} πειράματα στον φάκελο\n{path}")

    def clear_layout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget() is not None:
                child.widget().deleteLater()
            elif child.layout() is not None:
                self.clear_layout(child.layout())


app = QtWidgets.QApplication(sys.argv)
w = MainWindow()
app.exec_()
