# coding: utf-8

'''
Minimal example for PyMOL + LEAP Motion

(c) 2013 Schrödinger, Inc.

License: BSD-2-Clause
'''

import sys
import math
from pymol import cmd

if sys.platform.startswith('linux'):
    sys.path.append('/usr/share/Leap/LeapSDK/lib')
    sys.path.append('/usr/share/Leap/LeapSDK/lib/x64')
else:
    # TODO: ADJUST THIS
    sys.path.append('/Users/novak/Desktop/LeapDeveloperKit_2.2.0+23485_mac/LeapSDK/lib')

import Leap
from Leap import Matrix, Vector

class PymolListener(Leap.Listener):
    def __init__(self, target):
        super(PymolListener, self).__init__()

        self.target = target
        self.controller = Leap.Controller()
        self.controller.add_listener(self)


    def __del__(self):
        self.controller.remove_listener(self)

        super(PymolListener, self).__del__()

    def on_frame(self, controller):
        frame = controller.frame()
        self.target.send(frame)

def coroutine(func):
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start

@coroutine
def update_view():
    previous_frame = None
    while True:
        frame = (yield)
        if previous_frame:
            view = list(cmd.get_view())

            if frame.rotation_probability(previous_frame) > 0.1:
                m = frame.rotation_matrix(previous_frame)
                m *= Matrix(Vector(*view[0:3]),
                            Vector(*view[3:6]),
                            Vector(*view[6:9]))
                view[:9] = m.to_array_3x3()

            if frame.scale_probability(previous_frame) > 0.1:
                s = frame.scale_factor(previous_frame)
                delta_z = math.log(s) * 100.0
                view[11] += delta_z
                view[15] -= delta_z
                view[16] -= delta_z

            cmd.set_view(view)
        previous_frame = frame


listener = PymolListener(update_view())
