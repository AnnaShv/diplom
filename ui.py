# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui.ui'
#
# Created by: PyQt5 UI code generator 5.14.2
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1075, 280)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.groupBox_big1 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_big1.setObjectName("groupBox_big1")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.groupBox_big1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.groupBox_medium1 = QtWidgets.QGroupBox(self.groupBox_big1)
        self.groupBox_medium1.setTitle("")
        self.groupBox_medium1.setObjectName("groupBox_medium1")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.groupBox_medium1)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.groupBox_small11 = QtWidgets.QGroupBox(self.groupBox_medium1)
        self.groupBox_small11.setTitle("")
        self.groupBox_small11.setObjectName("groupBox_small11")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_small11)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_1 = QtWidgets.QLabel(self.groupBox_small11)
        self.label_1.setObjectName("label_1")
        self.verticalLayout_2.addWidget(self.label_1)
        self.label_2 = QtWidgets.QLabel(self.groupBox_small11)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.label_3 = QtWidgets.QLabel(self.groupBox_small11)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_2.addWidget(self.label_3)
        self.horizontalLayout_2.addWidget(self.groupBox_small11)
        self.groupBox_small12 = QtWidgets.QGroupBox(self.groupBox_medium1)
        self.groupBox_small12.setTitle("")
        self.groupBox_small12.setObjectName("groupBox_small12")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_small12)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.lineEdit_1 = QtWidgets.QLineEdit(self.groupBox_small12)
        self.lineEdit_1.setText("")
        self.lineEdit_1.setObjectName("lineEdit_1")
        self.verticalLayout_3.addWidget(self.lineEdit_1)
        self.lineEdit_2 = QtWidgets.QLineEdit(self.groupBox_small12)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.verticalLayout_3.addWidget(self.lineEdit_2)
        self.lineEdit_3 = QtWidgets.QLineEdit(self.groupBox_small12)
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.verticalLayout_3.addWidget(self.lineEdit_3)
        self.horizontalLayout_2.addWidget(self.groupBox_small12)
        self.horizontalLayout.addWidget(self.groupBox_medium1)
        self.groupBox_medium2 = QtWidgets.QGroupBox(self.groupBox_big1)
        self.groupBox_medium2.setTitle("")
        self.groupBox_medium2.setObjectName("groupBox_medium2")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.groupBox_medium2)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.groupBox_small21 = QtWidgets.QGroupBox(self.groupBox_medium2)
        self.groupBox_small21.setTitle("")
        self.groupBox_small21.setObjectName("groupBox_small21")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox_small21)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.label_4 = QtWidgets.QLabel(self.groupBox_small21)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_4.addWidget(self.label_4)
        self.label_5 = QtWidgets.QLabel(self.groupBox_small21)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_4.addWidget(self.label_5)
        self.label = QtWidgets.QLabel(self.groupBox_small21)
        self.label.setObjectName("label")
        self.verticalLayout_4.addWidget(self.label)
        self.horizontalLayout_3.addWidget(self.groupBox_small21)
        self.groupBox_small22 = QtWidgets.QGroupBox(self.groupBox_medium2)
        self.groupBox_small22.setTitle("")
        self.groupBox_small22.setObjectName("groupBox_small22")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox_small22)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.lineEdit_4 = QtWidgets.QLineEdit(self.groupBox_small22)
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.verticalLayout_5.addWidget(self.lineEdit_4)
        self.lineEdit_5 = QtWidgets.QLineEdit(self.groupBox_small22)
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.verticalLayout_5.addWidget(self.lineEdit_5)
        self.pushButton = QtWidgets.QPushButton(self.groupBox_small22)
        self.pushButton.setObjectName("pushButton")
        self.verticalLayout_5.addWidget(self.pushButton)
        self.horizontalLayout_3.addWidget(self.groupBox_small22)
        self.horizontalLayout.addWidget(self.groupBox_medium2)
        self.gridLayout.addWidget(self.groupBox_big1, 1, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.groupBox_big1.setTitle(_translate("MainWindow", "Дані"))
        self.label_1.setText(_translate("MainWindow", "Введіть кількість мешканців. "))
        self.label_2.setText(_translate("MainWindow", "Яку кількість днів хворий перебував поза ізолятором?"))
        self.label_3.setText(_translate("MainWindow", "Яку кількість днів хворий перебуває в ізоляторі?"))
        self.label_4.setText(_translate("MainWindow", "Інкубаційний період коронавірусу."))
        self.label_5.setText(_translate("MainWindow", "Індекс захворюваності."))
        self.label.setText(_translate("MainWindow", "Виконати"))
        self.pushButton.setText(_translate("MainWindow", "Натисни!"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())