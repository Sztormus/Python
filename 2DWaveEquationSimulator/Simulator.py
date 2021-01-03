import matplotlib
import numpy as np
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import QRegExp
from PyQt5.QtGui import QIntValidator, QRegExpValidator
from PyQt5.QtWidgets import QWidget, QLineEdit, QFormLayout, QLabel, QComboBox, QMessageBox, QPushButton
from Canvas import Canvas
from SimulationData import SimulationData

matplotlib.use('Qt5Agg')


class Simulator(QWidget):
    def __init__(self, *args, **kwargs):
        super(Simulator, self).__init__(*args, **kwargs)
        self.__simulation = SimulationData()

        self.__manage_window()

        self.__create_layouts()

        self.__create_buttons_and_inputs()

        self.__create_canvas()

        self.__manage_layouts()

        self.__connect_buttons()

        self.__set_timer()

        self.__change_initial_configuration()

        self.show()

    def __manage_window(self):
        self.setWindowTitle("2D Wave Equation Simulator")
        self.setFixedSize(1000, 860)

    def __create_layouts(self):
        self.__main_layout = QtWidgets.QVBoxLayout(self)
        self.__form_layout = QFormLayout(self)
        self.__start_stop_layout = QFormLayout(self)

    def __create_buttons_and_inputs(self):
        self.__start_button = QtWidgets.QPushButton("Start simulation", self)
        self.__stop_button = QtWidgets.QPushButton("Stop simulation", self)
        self.__start_button.setFixedWidth(485)
        self.__stop_button.setFixedWidth(485)
        self.__stop_button.setEnabled(False)

        self.__only_float = QRegExpValidator(QRegExp("-?[0-9]+\.?[0-9]+"))
        self.__only_float_time = QRegExpValidator(QRegExp("[0-9]+\.?[0-9]+"))
        self.__only_int = QIntValidator()
        self.__only_int.setBottom(0)
        self.__only_int.setTop(99)

        self.__label1 = QLabel("Number of points (X-axis):")
        self.__line1 = QLineEdit()
        self.__line1.setValidator(self.__only_int)
        self.__line1.setText("25")

        self.__label2 = QLabel("Number of points (Y-axis):")
        self.__line2 = QLineEdit()
        self.__line2.setValidator(self.__only_int)
        self.__line2.setText("25")

        self.__label3 = QLabel("Maximum simulation time [s]:")
        self.__line3 = QLineEdit()
        self.__line3.setValidator(self.__only_float_time)
        self.__line3.setText("1.0")

        self.__label4 = QPushButton("Alpha:")
        self.__label4.setFixedWidth(135)
        self.__line4 = QLineEdit()
        self.__line4.setValidator(self.__only_float)
        self.__line4.setText("0.0")

        self.__label5 = QPushButton("Beta:")
        self.__label5.setFixedWidth(135)
        self.__line5 = QLineEdit()
        self.__line5.setValidator(self.__only_float)
        self.__line5.setText("0.25")

        self.__label6 = QPushButton("Gamma:")
        self.__label6.setFixedWidth(135)
        self.__line6 = QLineEdit()
        self.__line6.setValidator(self.__only_float)
        self.__line6.setText("0.5")

        self.__label7 = QLabel("Fixed boundary conditions:")
        self.__line7 = QComboBox()
        self.__line7.addItems(["True", "False"])

        self.__label8 = QLabel("Zero-displacement initial condition:")
        self.__line8 = QComboBox()
        self.__line8.addItems(["False", "True"])

        self.__label9 = QLabel("Forcing term:")
        self.__line9 = QComboBox()
        self.__line9.addItems(["False", "True"])

    def __create_canvas(self):
        self.__canvas = Canvas(self)
        self.__canvas.setFixedSize(980, 566)

    def __manage_layouts(self):
        self.__form_layout.addRow(self.__label1, self.__line1)
        self.__form_layout.addRow(self.__label2, self.__line2)
        self.__form_layout.addRow(self.__label3, self.__line3)
        self.__form_layout.addRow(self.__label4, self.__line4)
        self.__form_layout.addRow(self.__label5, self.__line5)
        self.__form_layout.addRow(self.__label6, self.__line6)
        self.__form_layout.addRow(self.__label7, self.__line7)
        self.__form_layout.addRow(self.__label8, self.__line8)
        self.__form_layout.addRow(self.__label9, self.__line9)

        self.__start_stop_layout.addRow(self.__start_button, self.__stop_button)

        self.__main_layout.addWidget(self.__canvas)
        self.__main_layout.addLayout(self.__form_layout)
        self.__main_layout.addLayout(self.__start_stop_layout)

    def __connect_buttons(self):
        self.__start_button.clicked.connect(lambda: self.__start_simulation())
        self.__stop_button.clicked.connect(lambda: self.__stop_simulation())
        self.__label4.clicked.connect(lambda: self.__show_alpha_description())
        self.__label5.clicked.connect(lambda: self.__show_beta_description())
        self.__label6.clicked.connect(lambda: self.__show_gamma_description())
        self.__line1.textChanged.connect(lambda: self.__change_initial_configuration())
        self.__line2.textChanged.connect(lambda: self.__change_initial_configuration())
        self.__line3.textChanged.connect(lambda: self.__change_initial_configuration())
        self.__line4.textChanged.connect(lambda: self.__change_initial_configuration())
        self.__line5.textChanged.connect(lambda: self.__change_initial_configuration())
        self.__line6.textChanged.connect(lambda: self.__change_initial_configuration())
        self.__line7.currentTextChanged.connect(lambda: self.__change_initial_configuration())
        self.__line8.currentTextChanged.connect(lambda: self.__change_initial_configuration())
        self.__line9.currentTextChanged.connect(lambda: self.__change_initial_configuration())

    def __show_alpha_description(self):
        title = "Description of Alpha parameter"
        description = "Alpha is the factor responsible for damping the membrane movement"
        self.__show_message(QMessageBox.Information, title, description)

    def __show_beta_description(self):
        title = "Description of Beta parameter"
        description = "Beta is the scalar parameter of the Newmark algorithm"
        self.__show_message(QMessageBox.Information, title, description)

    def __show_gamma_description(self):
        title = "Description of Gamma parameter"
        description = "Gamma is the scalar parameter of the Newmark algorithm"
        self.__show_message(QMessageBox.Information, title, description)

    def __show_message(self, message_type, title, description):
        self.__msg = QMessageBox()
        self.__msg.setIcon(message_type)
        self.__msg.setWindowTitle(title)
        self.__msg.setText(description)
        self.__msg.exec_()

    def __disable_gui(self):
        self.__stop_button.setEnabled(True)
        self.__start_button.setEnabled(False)
        self.__line1.setEnabled(False)
        self.__line2.setEnabled(False)
        self.__line3.setEnabled(False)
        self.__line4.setEnabled(False)
        self.__line5.setEnabled(False)
        self.__line6.setEnabled(False)
        self.__line7.setEnabled(False)
        self.__line8.setEnabled(False)
        self.__line9.setEnabled(False)
        self.__label4.setEnabled(False)
        self.__label5.setEnabled(False)
        self.__label6.setEnabled(False)

    def __enable_gui(self):
        self.__stop_button.setEnabled(False)
        self.__start_button.setEnabled(True)
        self.__line1.setEnabled(True)
        self.__line2.setEnabled(True)
        self.__line3.setEnabled(True)
        self.__line4.setEnabled(True)
        self.__line5.setEnabled(True)
        self.__line6.setEnabled(True)
        self.__line7.setEnabled(True)
        self.__line8.setEnabled(True)
        self.__line9.setEnabled(True)
        self.__label4.setEnabled(True)
        self.__label5.setEnabled(True)
        self.__label6.setEnabled(True)

    def __set_timer(self):
        self.__timer = QtCore.QTimer()
        self.__timer.setInterval(100)
        self.__timer.timeout.connect(self.__continue_simulation)

    def __change_initial_configuration(self):
        try:
            self.__n_x = int(self.__line1.text())
            self.__n_y = int(self.__line2.text())
        except:
            return

        if self.__n_x < 3 or self.__n_y < 3:
            return

        if self.__line8.currentText() == "False":
            self.__zero_diaphragm_displacement = False
        elif self.__line8.currentText() == "True":
            self.__zero_diaphragm_displacement = True

        self.__x_net = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)
        self.__y_net = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)

        for j in range(self.__n_y):
            for i in range(self.__n_x):
                self.__x_net[j * self.__n_x + i] = i * 1.0 / self.__n_x
                self.__y_net[j * self.__n_x + i] = j * 1.0 / self.__n_y

        first_step = self.__simulation.get_first_step(self.__n_x, self.__n_y, self.__zero_diaphragm_displacement)

        self.__max_value = first_step.max() + 0.0001
        self.__min_value = -self.__max_value

        self.__canvas.axes.cla()
        self.__canvas.axes.set_xlabel("x [m]")
        self.__canvas.axes.set_ylabel("y [m]")
        self.__canvas.axes.set_zlabel("u [m]")
        if not self.__zero_diaphragm_displacement:
            self.__canvas.axes.set_zlim(self.__min_value, self.__max_value)
        else:
            self.__canvas.axes.set_zlim(-0.075, 0.075)

        self.__canvas.axes.set_title("Diaphragm displacement diagram u(x, y, t) \n t = 0.0 s")
        self.__canvas.axes.scatter(self.__x_net, self.__y_net, first_step, linewidth=0,
                                   antialiased=False)
        self.__canvas.draw()

    def __start_simulation(self):
        self.__n_x = int(self.__line1.text())
        self.__n_y = int(self.__line2.text())

        if self.__n_x < 3 or self.__n_y < 3:
            title = "Warning"
            description = "Number of points for X-axis and Y-axis must be at least 3"
            self.__show_message(QMessageBox.Warning, title, description)
            return

        self.__t = float(self.__line3.text())
        self.__alpha = float(self.__line4.text())
        self.__beta = float(self.__line5.text())
        self.__gamma = float(self.__line6.text())

        if self.__line7.currentText() == "False":
            self.__fixed_boundary_conditions = False
        elif self.__line7.currentText() == "True":
            self.__fixed_boundary_conditions = True

        if self.__line8.currentText() == "False":
            self.__zero_diaphragm_displacement = False
        elif self.__line8.currentText() == "True":
            self.__zero_diaphragm_displacement = True

        if self.__line9.currentText() == "False":
            self.__forcing_term_check = False
        elif self.__line9.currentText() == "True":
            self.__forcing_term_check = True

        self.__delta_t = 0.01
        self.__n_t = int(self.__t / self.__delta_t)
        self.__counter = 1

        self.__simulation.change_initial_configuration(self.__n_x, self.__n_y, self.__t, self.__alpha, self.__beta,
                                                       self.__gamma, self.__fixed_boundary_conditions,
                                                       self.__zero_diaphragm_displacement, self.__forcing_term_check)

        self.__disable_gui()

        self.__simulation.start_simulation()

        self.__canvas.axes.cla()
        self.__canvas.axes.set_xlabel("x [m]")
        self.__canvas.axes.set_ylabel("y [m]")
        self.__canvas.axes.set_zlabel("u [m]")
        if not self.__zero_diaphragm_displacement:
            self.__canvas.axes.set_zlim(self.__min_value, self.__max_value)
        else:
            self.__canvas.axes.set_zlim(-0.075, 0.075)

        self.__canvas.axes.set_title("Diaphragm displacement diagram u(x, y, t) \n t = 0.0 s")
        self.__canvas.axes.scatter(self.__x_net, self.__y_net, self.__simulation.get_actual_step(), linewidth=0,
                                   antialiased=False)
        self.__canvas.draw()

        self.__timer.start()

    def __continue_simulation(self):
        self.__canvas.axes.cla()
        self.__canvas.axes.set_xlabel("x [m]")
        self.__canvas.axes.set_ylabel("y [m]")
        self.__canvas.axes.set_zlabel("u [m]")
        if not self.__forcing_term_check:
            self.__canvas.axes.set_zlim(self.__min_value, self.__max_value)
        else:
            if self.__fixed_boundary_conditions or self.__simulation.get_actual_step().mean() < 0.05:
                self.__canvas.axes.set_zlim(-0.075, 0.075)
        self.__canvas.axes.set_title("Diaphragm displacement diagram u(x, y, t) \n t = "
                                     + str(round(self.__counter * self.__delta_t, 2)) + " s")
        self.__canvas.axes.scatter(self.__x_net, self.__y_net, self.__simulation.get_actual_step(), linewidth=0,
                                   antialiased=False)
        self.__canvas.draw()
        if self.__counter < self.__n_t:
            self.__counter += 1
            self.__simulation.calculate_next_step(self.__counter)
        else:
            self.__stop_simulation()

    def __stop_simulation(self):
        self.__timer.stop()
        self.__enable_gui()
