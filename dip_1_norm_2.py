import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy
import ui
from scipy.integrate import *
from pylab import *
import matplotlib.patches as mpatches
'''class PlotCanvas(FigureCanvas):
    def __init__(self, result, parent = None, width =5, height=4, dpi=100):
        fig=Figure(figsize=(width,height),dpi=dpi)
        self.axes=fig.add_subplot(111)
        FigureCanvas.__init__(self,fig)
        #self.setParent(parent)
        FigureCanvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot(result)
    def plot(self,result):
        data = result
        ax=self.figure.add_subplot(111)
        ax.plot(data,'r-')
        ax.set-title('Z')
        self.draw()'''


#$blabla 000


class App(QtWidgets.QMainWindow, ui.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        title = 'Прогнозування захворювання COVID-19'
        self.setWindowTitle(title)

        '''попытка
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)'''
        self.pushButton.clicked.connect(self.calc)

        #m = PlotCanvas(self)
        #self.show()
        #self.pushButton.clicked.connect(self.execute)
        #


    def plot(self, result):
        '''data = [result[:,2]]
        ax = self.figure.add_subplot(111)
        ax.clear()
        ax.plot(data, '*-')
        self.canvas.draw()'''
        #plot(t,S,color='green',lw=2)
        xlabel('Кількість днів')
        ylabel("Кількість людей")

        # plot(t,E,col"or='blue',lw=2)
        t = linspace(0, 100, 100)
        S = result[:, 0]
        E = result[:, 1]
        Z = result[:, 2]
        I = result[:, 3]
        R = result[:, 4]
        D = result[:, 5]
        ZI = I + Z

        plot(t, Z, color='red', lw=2)

        plot(t, I, color='green', lw=2)

        # plot(t,R,color='magenta',lw=2)

        plot(t, D, color='grey', lw=2)
        plot(t, ZI, color='black', lw=2)

        xlim(0)

        # patch_S = mpatches.Patch(color='green', label='Сприйнятливі')
        # patch_E = mpatches.Patch(color='blue', label='Латентні')
        patch_Z = mpatches.Patch(color='red', label='Заразні')
        patch_I = mpatches.Patch(color='green', label='Ізольовані хворі')
        # patch_R = mpatches.Patch(color='magenta', label='Одужавші')
        patch_D = mpatches.Patch(color='grey', label='Померлі')
        patch_ZI = mpatches.Patch(color='black', label='Всі хворі')
        plt.legend(handles=[patch_Z, patch_I, patch_D, patch_ZI])
        grid()
        show()


    def fun777(self, y, t, R0, T2, N, T1, T3):
            S, E, Z, I, R, D = y
            return [(-R0 / (T2 * N)) * S * Z,
                    (R0 / (T2 * N)) * S * Z - 1 / T1 * E,
                    1 / T1 * E - 1 / T2 * Z,
                    1 / T2 * Z - 1 / T3 * I,
                    (1 - 0.001) / T3 * Z,
                    0.001 / T3 * Z]

    def calc(self):
        t = linspace(0, 100, 100)
        '''R0=30
        T2=5
        T1=5.2
        T3=5
        N=50
        z=1'''


        N = int(self.lineEdit_1.text())
        S0, E0, Z0, I0, R0, D0 = N - 1, 0, 1, 0, 0, 0
        y0 = S0, E0, Z0, I0, R0, D0

        '''def f0(y,t):
            S, E, Z, R, D = y
            S, E, Z, R, D = 999,0,1,0,0
            t=0
            return [-1.25*S(t)*Z(t), 1.25*S(t)*Z(t)-1/5.2*E(t), 1/5.2*E(t)-0.1*Z(t), 0.0999*Z(t), 0.0001*Z(t)]'''
        R0 = int(self.lineEdit_5.text())
        T2 = float(self.lineEdit_2.text())

        T1 = float(self.lineEdit_4.text())
        T3 = float(self.lineEdit_3.text())
        result = odeint(self.fun777, y0, t, args=(R0, T2, N, T1, T3,))
        S = result[:, 0]
        E = result[:, 1]
        Z = result[:, 2]
        I = result[:, 3]
        R = result[:, 4]
        D = result[:, 5]
        ZI = I + Z
        # plot(t,S,color='green',lw=2)
        xlabel('Кількість днів')
        ylabel("Кількість людей")

        # plot(t,E,color='blue',lw=2)

        plot(t, Z, color='red', lw=2)

        plot(t, I, color='green', lw=2)

        # plot(t,R,color='magenta',lw=2)

        plot(t, D, color='grey', lw=2)
        plot(t, ZI, color='black', lw=2)

        xlim(0)

        # patch_S = mpatches.Patch(color='green', label='Сприйнятливі')
        # patch_E = mpatches.Patch(color='blue', label='Латентні')
        patch_Z = mpatches.Patch(color='red', label='Заразні')
        patch_I = mpatches.Patch(color='green', label='Ізольовані хворі')
        # patch_R = mpatches.Patch(color='magenta', label='Одужавші')
        patch_D = mpatches.Patch(color='grey', label='Померлі')
        patch_ZI = mpatches.Patch(color='black', label='Всі хворі')
        plt.legend(handles=[patch_Z, patch_I, patch_D, patch_ZI])
        grid()
        show()
        #self.plot(result)
        #m = PlotCanvas(self,Z, width=5, height=4)
        #m.move(0,0)
        #self.show()

def main():
    app=QtWidgets.QApplication(sys.argv)
    window=App()
    window.show()
    app.exec_()

if __name__=="__main__":
    main()

#kva kva

'''
print("Введите количество жителей общежития")
N = input()
N=int(N)
print("Введите количество зараженных жителей")
z = input()
z=int(z)
print("Индекс заражаемости COVID-19")
R0 = input()
R0=float(R0)
print("Сколько времени больной не придерживался карантина?")
T2 = input()
T2 = int(T2)
print("Инкубационный период COVID-19")
T1 = input()
T1 = float(T1)'''

'''print("Количество жителей: {0}".format(N))'''
'''
t = linspace(0,100,100)
#R0=30
#T2=5
#T1=5.2
#T3=5
#N=50
#z=1
def f(y,t, R0,T2,N,T1,T3):
    S,E,Z,I,R,D=y
    return [(-R0/(T2*N))*S*Z,
            (R0/(T2*N))*S*Z-1/T1*E,
            1/T1*E-1/T2*Z,
            1/T2*Z - 1/T3*I,
            (1-0.001)/T3*Z,
            0.001/T3*Z]
N=int(input("Населення гуртожитка/студмістечка? "))
S0,E0,Z0,I0,R0,D0=N-1,0,1,0,0,0
y0=S0,E0,Z0,I0,R0,D0

def f0(y,t):
    S, E, Z, R, D = y
    S, E, Z, R, D = 999,0,1,0,0
    t=0
    return [-1.25*S(t)*Z(t), 1.25*S(t)*Z(t)-1/5.2*E(t), 1/5.2*E(t)-0.1*Z(t), 0.0999*Z(t), 0.0001*Z(t)]
R0=int(input("Індекс захворюваності = "))
T2=float(input("Скільки днів хвора людина перебувала поза ізолятором? "))

T1=float(input("Інкубаційний період коронавірусу = "))
T3=float(input("Період, який людина перебуває в ізоляторі "))
result = odeint(f,y0,t,args=(R0,T2,N,T1,T3))
S = result[:,0]
E = result[:,1]
Z = result[:,2]
I = result[:,3]
R = result[:,4]
D = result[:,5]
ZI=I+Z'''

'''
#plot(t,S,color='green',lw=2)
xlabel('Кількість днів')
ylabel("Кількість людей")

#plot(t,E,color='blue',lw=2)

plot(t,Z,color='red',lw=2)

plot(t,I,color='green',lw=2)

#plot(t,R,color='magenta',lw=2)

plot(t,D,color='grey',lw=2)
plot(t,ZI,color='black',lw=2)

xlim(0)

#patch_S = mpatches.Patch(color='green', label='Сприйнятливі')
#patch_E = mpatches.Patch(color='blue', label='Латентні')
patch_Z = mpatches.Patch(color='red', label='Заразні')
patch_I = mpatches.Patch(color='green', label='Ізольовані хворі')
#patch_R = mpatches.Patch(color='magenta', label='Одужавші')
patch_D = mpatches.Patch(color='grey', label='Померлі')
patch_ZI = mpatches.Patch(color='black', label='Всі хворі')
plt.legend(handles=[ patch_Z,patch_I, patch_D, patch_ZI])
grid()
show()

def odeint():
         #dy1/dt=y2
         #dy2/dt=y1**2+1:
         def f(y,t):
                   return y**2+1
         t = arange(0,1,0.01)
         y0 =0.0
         y=odeint(f, y0,t)
         y = array(y).flatten()
         return y,t

def f(t, y):
         f = zeros([1])
         f[0] = y[0]**2+1
         return f
to = 0.
tEnd = 1
yo = array([0.])
tau = 0.01
z=1
y,t=odeint()
plt.plot(t,abs(array(tan(t))-array(y)),label='Функция odein')

plt.show()'''
