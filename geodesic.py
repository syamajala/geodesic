import os
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt'
import numpy as np
import sympy as sp
import networkx as nx
import itertools as it
from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
    SceneEditor

pi = np.pi
N = 50
eps = 1/N
BLACK = (0, 0, 0)


class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True)

    def __init__(self, surf, **traits):
        super(HasTraits, self).__init__(**traits)

        self.surf = surf
        self.plane = surf.plane

        self.f = sp.lambdify(self.surf.pvars, self.surf.parameterization,
                             [{'ImmutableMatrix': np.array}, "numpy"])

        self.figure = self.scene.mlab.gcf()
        self.figure2 = self.scene.mlab.figure(2)

        self.plane_points = []
        self.surf_points = []

        self.j = self.surf.parameterization.jacobian(self.pvars)
        self.G = (self.j.transpose()*self.j)
        self.G.simplify()
        self.Ginv = self.G.inv()

        # self.imesh = InducedMesh(10, self.f)

    def christoffel(self, k):
        return sp.Matrix([[self.christoffel_ijk(0, 0, k), self.christoffel_ijk(0, 1, k)],
                          [self.christoffel_ijk(1, 0, k), self.christoffel_ijk(1, 1, k)]])

    def christoffel_ijk(self, i, j, k):
        k = k - 1
        sum = 0

        for l in range(0, 2):
            sum += self.Ginv[k, l] * (sp.diff(self.G[j, l], self.surf.pvars[i]) +
                                      sp.diff(self.G[l, i], self.surf.pvars[j]) -
                                      sp.diff(self.G[i, j], self.surf.pvars[l]))

        return sum/2

    def curve(self, l):
        """Given a plane curve l, find the surface curve"""
        return [self.f(i[0], i[1]) for i in l]

    def y_pk(self, p, k, x):
        k = k - 1

        if p == 0 or p == N:
            return x[p][k]

        sum = 0
        for e in it.product(range(0, 2), range(0, 2)):
            i, j = e
            ch = self.christoffel_ijk(i, j, k)
            sum += ch.subs(zip(self.surf.pvars, x[p])) * \
                   (x[p+1][i] - x[p-1][i]) * (x[p+1][j] - x[p-1][j])

        sum = (x[p+1][k] + x[p-1][k])/2 + sum/4

    @on_trait_change('scene.activated')
    def update_plot(self):
        x, y, z = self.f(self.plane.x, self.plane.y)[0]
        self.scene.mlab.mesh(x, y, z, color=(0, 1, 0))

        self.scene.mlab.figure(self.figure2)
        self.scene.mlab.clf()

        picker = self.figure2.on_mouse_pick(self.onpick)
        picker.tolerance = 0.01

        self.scene.mlab.surf(self.plane.x, self.plane.y, self.plane.f())

    def onpick(self, event):
        point = event.pick_position
        added = False
        if self.plane.in_plane(point) and point not in self.plane_points:
            print "adding point",
            self.plane_points.append(point)
            self.scene.mlab.points3d(*point, scale_factor=0.1)
            self.scene.mlab.text3d(*point, scale=0.2, color=BLACK,
                                   text=str(len(self.plane_points)))

            added = True
        else:
            print "plane already contains point",

        print point

        l = None

        if len(self.plane_points) >= 2 and added:
            start = self.plane_points[-2]
            l = self.plane.line(start, point)
            x = []
            y = []
            z = []
            for i in l:
                x.append(i[0])
                y.append(i[1])
                z.append(i[2])

            self.scene.mlab.plot3d(x, y, z, tube_radius=0.01)

        if added:
            self.scene.mlab.figure(self.figure)
            spoint = self.surf.f(point[0], point[1])[0]
            self.surf_points.append(spoint)
            self.scene.mlab.points3d(*spoint, scale_factor=0.1)
            self.scene.mlab.text3d(*spoint, scale=0.2, color=BLACK,
                                   text=str(len(self.surf_points)))

            if l is not None:
                c = self.curve(l)

                x = []
                y = []
                z = []
                for i in c:
                    j = i[0]
                    x.append(j[0])
                    y.append(j[1])
                    z.append(j[2])

                self.scene.mlab.plot3d(x, y, z, tube_radius=0.01, color=BLACK)

            self.scene.mlab.figure(self.figure2)


class Plane():

    def __init__(self, x1, xs, y1, ys, x0=0, y0=0, z=0):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.z = z
        self.x, self.y = np.mgrid[x0:x1:xs, y0:y1:ys]

    def f(self):
        return np.array([[self.z]*self.x.shape[1]]*self.x.shape[0])

    def in_plane(self, p):
        """Checks if the point p is contained in the plane."""
        x, y, z = p
        if self.x0 <= x and x <= self.x1:
            if self.y0 <= y and y <= self.y1:
                if self.z == z:
                    return True

        return False

    def line(self, p, q):
        """Line from p to q, ie. p + t*(q-p)"""
        t = np.linspace(0, 1, N)
        p_ = np.array(p)
        q_ = np.array(q)

        return np.array([p_ + i*(q_-p_) for i in t])


class InducedMesh():

    def __init__(self, size, f):
        self.size = size
        self.f = f

        self.x, self.y = np.mgrid[0:1:1j*self.size, 0:1:1j*self.size]
        self.w = self.f(self.x, self.y)[0]

        self.graph = nx.Graph()
        self.graph.add_nodes_from(range(0, self.size**2))


class Sphere(Visualization):

    def __init__(self, radius):

        self.plane = Plane(pi, 101j, 2 * pi, 101j)
        self.radius = radius

        self.theta, self.phi = sp.symbols("theta phi")
        self.pvars = sp.Matrix([self.theta, self.phi])
        self.parameterization = sp.Matrix([[self.radius * sp.sin(self.phi) *
                                            sp.cos(self.theta),
                                            self.radius * sp.sin(self.phi) *
                                            sp.sin(self.theta),
                                            self.radius * sp.cos(self.phi)]])

        super(Sphere, self).__init__(self)

    # @on_trait_change('scene.activated')
    # def update_plot(self):
    #     super(Sphere, self).update_plot()


class Torus(Visualization):

    def __init__(self, R, r):

        self.plane = Plane(2 * pi, 101j, 2 * pi, 101j)
        self.r = r
        self.R = R

        self.theta, self.phi = sp.symbols("theta phi")
        self.pvars = sp.Matrix([self.theta, self.phi])
        self.parameterization = sp.Matrix(
            [[(self.R + self.r*sp.cos(self.phi))*sp.cos(self.theta),
              (self.R + self.r*sp.cos(self.phi))*sp.sin(self.theta),
              self.r*sp.sin(self.phi)]])

        super(Torus, self).__init__(self)


class MayaviQWidget(QtGui.QWidget):
    def __init__(self, visualization, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = visualization
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)


if __name__ == "__main__":
    app = QtGui.QApplication.instance()
    container = QtGui.QWidget()
    container.setWindowTitle("Embedding Mayavi in a PyQt4 Application")
    layout = QtGui.QGridLayout(container)
#    s = Sphere(1)
    s = Torus(3, 2)
    mayavi_widget = MayaviQWidget(s, container)
    layout.addWidget(mayavi_widget, 1, 1)
    label = QtGui.QLabel(container)
    label.setText("hi")
    layout.addWidget(label, 1, 2)
    container.show()
    window = QtGui.QMainWindow()
    window.setCentralWidget(container)
    window.show()
    app.exec_()
else:
    s = Sphere(1)
    print s.christoffel(1)
