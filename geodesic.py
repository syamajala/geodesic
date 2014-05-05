import os
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt'
import numpy as np
import sympy as sp
from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
    SceneEditor

pi = np.pi


class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True)

    def __init__(self, surf, **traits):
        super(HasTraits, self).__init__(**traits)

        self.surf = surf
        self.plane = surf.plane
        self.figure = self.scene.mlab.gcf()
        self.figure2 = self.scene.mlab.figure(2)

        self.plane_points = []
        self.surf_points = []

    def plane(self, x, y):
        return np.array([[0]*x.shape[1]]*x.shape[0])

    @on_trait_change('scene.activated')
    def update_plot(self):
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
            added = True
        else:
            print "plane already contains point",

        print point

        if len(self.plane_points) >= 2:
            start = self.plane_points[-2]
            x = [start[0], point[0]]
            y = [start[1], point[1]]
            z = [start[2], point[2]]
            self.scene.mlab.plot3d(x, y, z, tube_radius=0.01)

        if added == True:
            self.scene.mlab.figure(self.figure)
            self.scene.mlab.points3d(*self.surf.f(point[0], point[1])[0],
                                     scale_factor=0.1)


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
        x, y, z = p
        if self.x0 <= x and x <= self.x1:
            if self.y0 <= y and y <= self.y1:
                if self.z == z:
                    return True

        return False


class Sphere(Visualization):

    def __init__(self, radius):
        super(Sphere, self).__init__(self)
        self.plane = Plane(pi, 101j, 2 * pi, 101j)
        self.radius = radius

        self.phi, self.theta = sp.symbols("phi theta")
        self.sphere = sp.Matrix([[self.radius * sp.sin(self.phi) * sp.cos(self.theta),
                                  self.radius * sp.sin(self.phi) * sp.sin(self.theta),
                                  self.radius * sp.cos(self.phi)]])
        self.f = sp.lambdify((self.phi, self.theta), self.sphere,
                             [{'ImmutableMatrix': np.array}, "numpy"])

    @on_trait_change('scene.activated')
    def update_plot(self):
        x, y, z = self.f(self.plane.x, self.plane.y)[0]
        self.scene.mlab.mesh(x, y, z, color=(0, 1, 0))

        super(Sphere, self).update_plot()


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
    s = Sphere(1)
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
