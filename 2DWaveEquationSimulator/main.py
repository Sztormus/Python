import sys
from PyQt5.QtWidgets import QApplication
from Simulator import Simulator

app = QApplication(sys.argv)
w = Simulator()
app.exec_()
