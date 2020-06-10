# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from Params import Params
from sol import Solution
import matplotlib.pyplot as plt
import numpy as np
import sys


class Default_params:
    J = None
    Jz = None
    l0 = None
    m = None
    c = None
    s0 = None
    omega = None
    eps = None
    step = None
    alpha = None
    time = None

    def __init__(self, J=0.74, Jz=0.017, l0=1, m=13.7, c=20, s0=0.0005, omega=314, eps=0.1, step=0.1, alpha=0.5, time=10):
        self.J = J
        self.Jz = Jz
        self.l0 = l0
        self.m = m
        self.c = c
        self.s0 = s0
        self.omega = omega
        self.eps = eps
        self.step = step
        self.alpha = alpha
        self.time = time


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1109, 911)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(0, 0, 411, 561))
        self.groupBox.setObjectName("groupBox")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setGeometry(QtCore.QRect(60, 30, 31, 21))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setGeometry(QtCore.QRect(60, 70, 111, 21))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setGeometry(QtCore.QRect(310, 20, 91, 31))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.groupBox)
        self.label_5.setGeometry(QtCore.QRect(310, 60, 71, 31))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.groupBox)
        self.label_6.setGeometry(QtCore.QRect(50, 100, 31, 31))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.groupBox)
        self.label_7.setGeometry(QtCore.QRect(330, 100, 61, 21))
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(self.groupBox)
        self.label_8.setGeometry(QtCore.QRect(50, 140, 41, 31))
        self.label_8.setObjectName("label_8")
        self.label_9 = QtWidgets.QLabel(self.groupBox)
        self.label_9.setGeometry(QtCore.QRect(330, 140, 51, 16))
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.groupBox)
        self.label_10.setGeometry(QtCore.QRect(50, 180, 41, 31))
        self.label_10.setObjectName("label_10")
        self.label_11 = QtWidgets.QLabel(self.groupBox)
        self.label_11.setGeometry(QtCore.QRect(330, 180, 41, 16))
        self.label_11.setObjectName("label_11")
        self.label_12 = QtWidgets.QLabel(self.groupBox)
        self.label_12.setGeometry(QtCore.QRect(30, 220, 101, 31))
        self.label_12.setObjectName("label_12")
        self.label_13 = QtWidgets.QLabel(self.groupBox)
        self.label_13.setGeometry(QtCore.QRect(320, 220, 61, 21))
        self.label_13.setObjectName("label_13")
        self.label_14 = QtWidgets.QLabel(self.groupBox)
        self.label_14.setGeometry(QtCore.QRect(10, 280, 171, 16))
        self.label_14.setObjectName("label_14")
        self.label_15 = QtWidgets.QLabel(self.groupBox)
        self.label_15.setGeometry(QtCore.QRect(10, 320, 181, 21))
        self.label_15.setObjectName("label_15")
        self.label_18 = QtWidgets.QLabel(self.groupBox)
        self.label_18.setGeometry(QtCore.QRect(30, 360, 91, 31))
        self.label_18.setObjectName("label_18")
        self.a = QtWidgets.QLabel(self.groupBox)
        self.a.setGeometry(QtCore.QRect(50, 400, 51, 31))
        self.a.setObjectName("a")
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setGeometry(QtCore.QRect(50, 450, 61, 21))
        self.label.setObjectName("label")
        self.label_16 = QtWidgets.QLabel(self.groupBox)
        self.label_16.setGeometry(QtCore.QRect(320, 450, 61, 21))
        self.label_16.setObjectName("label_16")
        self.J = QtWidgets.QLineEdit(self.groupBox)
        self.J.setGeometry(QtCore.QRect(160, 20, 111, 31))
        self.J.setObjectName("J")
        self.Jz = QtWidgets.QLineEdit(self.groupBox)
        self.Jz.setGeometry(QtCore.QRect(160, 60, 113, 31))
        self.Jz.setObjectName("Jz")
        self.l0 = QtWidgets.QLineEdit(self.groupBox)
        self.l0.setGeometry(QtCore.QRect(160, 100, 111, 31))
        self.l0.setObjectName("l0")
        self.m = QtWidgets.QLineEdit(self.groupBox)
        self.m.setGeometry(QtCore.QRect(160, 140, 111, 31))
        self.m.setObjectName("m")
        self.s0 = QtWidgets.QLineEdit(self.groupBox)
        self.s0.setGeometry(QtCore.QRect(160, 180, 111, 31))
        self.s0.setObjectName("s0")
        self.omega = QtWidgets.QLineEdit(self.groupBox)
        self.omega.setGeometry(QtCore.QRect(160, 220, 113, 31))
        self.omega.setObjectName("omega")
        self.c = QtWidgets.QLineEdit(self.groupBox)
        self.c.setGeometry(QtCore.QRect(160, 270, 113, 31))
        self.c.setObjectName("c")
        self.eps = QtWidgets.QLineEdit(self.groupBox)
        self.eps.setGeometry(QtCore.QRect(160, 310, 113, 31))
        self.eps.setObjectName("eps")
        self.step = QtWidgets.QLineEdit(self.groupBox)
        self.step.setGeometry(QtCore.QRect(160, 360, 113, 31))
        self.step.setObjectName("step")
        self.alpha = QtWidgets.QLineEdit(self.groupBox)
        self.alpha.setGeometry(QtCore.QRect(160, 400, 113, 31))
        self.alpha.setObjectName("alpha")
        self.time = QtWidgets.QLineEdit(self.groupBox)
        self.time.setGeometry(QtCore.QRect(160, 440, 113, 31))
        self.time.setObjectName("time")
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setGeometry(QtCore.QRect(0, 580, 411, 281))
        self.groupBox_2.setObjectName("groupBox_2")
        self.radioButtonHarmonic = QtWidgets.QRadioButton(self.groupBox_2)
        self.radioButtonHarmonic.setGeometry(QtCore.QRect(10, 40, 161, 31))
        self.radioButtonHarmonic.setObjectName("radioButtonHarmonic")
        self.radioButtonDumped = QtWidgets.QRadioButton(self.groupBox_2)
        self.radioButtonDumped.setGeometry(QtCore.QRect(10, 80, 171, 31))
        self.radioButtonDumped.setObjectName("radioButtonDumped")
        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(550, 40, 301, 121))
        self.label_17.setObjectName("label_17")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        self.label_19.setGeometry(QtCore.QRect(550, 140, 231, 51))
        self.label_19.setObjectName("label_19")
        self.max_u = QtWidgets.QLabel(self.centralwidget)
        self.max_u.setGeometry(QtCore.QRect(800, 90, 221, 21))
        self.max_u.setText("")
        self.max_u.setObjectName("max_u")
        self.Max_l = QtWidgets.QLabel(self.centralwidget)
        self.Max_l.setGeometry(QtCore.QRect(800, 145, 231, 31))
        self.Max_l.setText("")
        self.Max_l.setObjectName("Max_l")
        self.label_20 = QtWidgets.QLabel(self.centralwidget)
        self.label_20.setGeometry(QtCore.QRect(550, 260, 61, 31))
        self.label_20.setObjectName("label_20")
        self.label_21 = QtWidgets.QLabel(self.centralwidget)
        self.label_21.setGeometry(QtCore.QRect(550, 310, 51, 31))
        self.label_21.setObjectName("label_21")
        self.label_22 = QtWidgets.QLabel(self.centralwidget)
        self.label_22.setGeometry(QtCore.QRect(550, 360, 51, 31))
        self.label_22.setObjectName("label_22")
        self.label_23 = QtWidgets.QLabel(self.centralwidget)
        self.label_23.setGeometry(QtCore.QRect(550, 410, 51, 31))
        self.label_23.setObjectName("label_23")
        self.maxu1 = QtWidgets.QLabel(self.centralwidget)
        self.maxu1.setGeometry(QtCore.QRect(610, 260, 221, 31))
        self.maxu1.setText("")
        self.maxu1.setObjectName("maxu1")
        self.maxu2 = QtWidgets.QLabel(self.centralwidget)
        self.maxu2.setGeometry(QtCore.QRect(610, 310, 221, 31))
        self.maxu2.setText("")
        self.maxu2.setObjectName("maxu2")
        self.maxu3 = QtWidgets.QLabel(self.centralwidget)
        self.maxu3.setGeometry(QtCore.QRect(610, 360, 221, 31))
        self.maxu3.setText("")
        self.maxu3.setObjectName("maxu3")
        self.maxu4 = QtWidgets.QLabel(self.centralwidget)
        self.maxu4.setGeometry(QtCore.QRect(610, 410, 221, 31))
        self.maxu4.setText("")
        self.maxu4.setObjectName("maxu4")
        self.pushButtonRun = QtWidgets.QPushButton(self.centralwidget)
        self.pushButtonRun.setGeometry(QtCore.QRect(550, 480, 181, 61))
        self.pushButtonRun.setObjectName("pushButtonRun")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1109, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.init_input_containers()
        self.pushButtonRun.clicked.connect(self.display_solution)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Rotor"))
        self.groupBox.setTitle(_translate("MainWindow", "Параметры"))
        self.label_2.setText(_translate("MainWindow", "J"))
        self.label_3.setText(_translate("MainWindow", "Jz"))
        self.label_4.setText(_translate("MainWindow", "кг *  м^2"))
        self.label_5.setText(_translate("MainWindow", "кг * м^2"))
        self.label_6.setText(_translate("MainWindow", "Длина"))
        self.label_7.setText(_translate("MainWindow", "м"))
        self.label_8.setText(_translate("MainWindow", "Масса"))
        self.label_9.setText(_translate("MainWindow", "кг"))
        self.label_10.setText(_translate("MainWindow", "s0"))
        self.label_11.setText(_translate("MainWindow", "м"))
        self.label_12.setText(_translate("MainWindow", "Угловая частота"))
        self.label_13.setText(_translate("MainWindow", "рад/с"))
        self.label_14.setText(_translate("MainWindow", "Коэффициент жесткости"))
        self.label_15.setText(_translate("MainWindow", "eps - Коэффициент трения"))
        self.label_18.setText(_translate("MainWindow", "Дискретный шаг"))
        self.a.setText(_translate("MainWindow", "alpha"))
        self.label.setText(_translate("MainWindow", "Время"))
        self.label_16.setText(_translate("MainWindow", "c"))
        self.groupBox_2.setTitle(_translate("MainWindow", "Внешнее воздействие"))
        self.radioButtonHarmonic.setText(_translate("MainWindow", "Гармонические колебания"))
        self.radioButtonDumped.setText(_translate("MainWindow", "Затухающие колебания"))
        self.label_17.setText(_translate("MainWindow", "Максимальное отклонение в верхней точке"))
        self.label_19.setText(_translate("MainWindow", "Максимальное отклонение в нижней точке"))
        self.label_20.setText(_translate("MainWindow", "Max u1"))
        self.label_21.setText(_translate("MainWindow", "Max u2"))
        self.label_22.setText(_translate("MainWindow", "Max u3"))
        self.label_23.setText(_translate("MainWindow", "Max u4"))
        self.pushButtonRun.setText(_translate("MainWindow", "Run"))

    def run(self):
        app = QtWidgets.QApplication(sys.argv)
        MainWindow = QtWidgets.QMainWindow()
        ui = Ui_MainWindow()
        ui.setupUi(MainWindow)
        MainWindow.show()
        sys.exit(app.exec_())

    def init_input_containers(self):
        default_params = Default_params()
        self.J.setText(str(default_params.J))
        self.Jz.setText(str(default_params.Jz))
        self.m.setText(str(default_params.m))
        self.l0.setText(str(default_params.l0))
        self.c.setText(str(default_params.c))
        self.omega.setText(str(default_params.omega))
        self.eps.setText(str(default_params.eps))
        self.step.setText(str(default_params.step))
        self.time.setText(str(default_params.time))
        self.s0.setText(str(default_params.s0))
        self.alpha.setText(str(default_params.alpha))
        self.radioButtonDumped.setChecked(True)

    def display_solution(self):
        end_time = float(self.time.text())
        J = float(self.J.text())
        Jz = float(self.Jz.text())
        m = float(self.m.text())
        l0 = float(self.l0.text()) / 4
        c = float(self.c.text())
        omega = float(self.omega.text())
        eps = float(self.eps.text())
        step = float(self.step.text())
        s0 = float(self.s0.text())
        alpha = float(self.alpha.text())
        is_harmonic = False

        t = np.arange(0, end_time + 0.001, 0.001)

        if self.radioButtonHarmonic.isChecked():
            is_harmonic = True

        params = Params(J=J, Jz=Jz, l0=l0, m=m, c=c, s0=s0, omega=omega, eps=eps, step=step, alpha=alpha)
        disc_t = np.arange(0, end_time + step, step)
        solution = Solution(end_time=end_time, params=params)
        z1_cont, z2_cont, z1_cont_not_norm = solution.get_continuous_solution(params=params, is_harmonic=is_harmonic)
        z1_disc, z2_disc, z1_disc_not_norm = solution.get_discrete_solution(params=params, is_harmonic=is_harmonic)
        z1, z2, z1_not_norm = solution.get_solution(params=params, is_harmonic=is_harmonic)

        self.max_u.setText(str(z1[0].max()))
        self.Max_l.setText(str(z1[1].max()))
        self.maxu1.setText(str(z2[0].max()))
        self.maxu2.setText(str(z2[1].max()))
        self.maxu3.setText(str(z2[2].max()))
        self.maxu4.setText(str(z2[3].max()))

        with open('log.txt', 'w') as f:
            f.write('z1_u: ')
            for k in range(len(z1[0])):
                f.write("{} ".format(z1[0][k]))
            f.write('\nz1_l: ')
            for k in range(len(z1[1])):
                f.write("{} ".format(z1[1][k]))

            f.write('\nu1: ')
            for k in range(len(z2[0])):
                f.write("{} ".format(z2[0][k]))
            f.write('\nu2: ')
            for k in range(len(z2[1])):
                f.write("{} ".format(z2[1][k]))
            f.write('\nu3: ')
            for k in range(len(z2[2])):
                f.write("{} ".format(z2[2][k]))
            f.write('\nu4: ')
            for k in range(len(z2[3])):
                f.write("{} ".format(z2[3][k]))

        plt.figure(num='Графики переходных процессов')
        plt.subplot(2, 3, 3)
        plt.plot(t, z1_cont[0], color='blue')
        plt.plot(disc_t, z1_disc[0], color='red', linestyle='--')
        plt.plot(t, z1[0], color='black', linestyle='-.')
        plt.ylabel('z1_u')
        plt.grid()
        plt.subplot(2, 3, 6)
        plt.plot(t, z1_cont[1], color='blue')
        plt.plot(disc_t, z1_disc[1], color='red', linestyle='--')
        plt.plot(t, z1[1], color='black', linestyle='-.')
        plt.ylabel('z1_l')
        plt.grid()
        plt.legend(['Continuous', 'Solution'])
        plt.subplot(2, 3, 1)
        plt.plot(t, z2_cont[0], color='blue')
        plt.plot(disc_t, z2_disc[0], color='red', linestyle='--')
        plt.plot(t, z2[0], color='black', linestyle='-.')
        plt.ylabel('u_1')
        plt.grid()
        plt.subplot(2, 3, 2)
        plt.plot(t, z2_cont[1], color='blue')
        plt.plot(disc_t, z2_disc[1], color='red', linestyle='--')
        plt.plot(t, z2[1], color='black', linestyle='-.')
        plt.ylabel('u_2')
        plt.grid()
        plt.subplot(2, 3, 4)
        plt.plot(t, z2_cont[2], color='blue')
        plt.plot(disc_t, z2_disc[2], color='red', linestyle='--')
        plt.plot(t, z2[2], color='black', linestyle='-.')
        plt.ylabel('u_3')
        plt.grid()
        plt.subplot(2, 3, 5)
        plt.plot(t, z2_cont[3], color='blue')
        plt.plot(disc_t, z2_disc[3], color='red', linestyle='--')
        plt.plot(t, z2[3], color='black', linestyle='-.')
        plt.ylabel('u_4')
        plt.grid()

        plt.figure(num='Смещения')
        plt.subplot(1, 2, 1)
        plt.plot(z1_not_norm[0], z1_not_norm[2])
        plt.title('Смещения в верхних подшипниках')
        plt.xlabel('x_u')
        plt.ylabel('y_u')
        plt.grid()
        plt.subplot(1, 2, 2)
        plt.plot(z1_not_norm[1], z1_not_norm[3])
        plt.title('Смещения в нижних подшипниках')
        plt.xlabel('x_l')
        plt.ylabel('y_l')
        plt.grid()
        plt.show()
