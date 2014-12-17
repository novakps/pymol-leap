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
    sys.path.append('/opt/LeapDeveloperKit_2.2.0+23485_mac/LeapSDK/lib')

import Leap
from Leap import Matrix, Vector

class PymolListener(Leap.Listener):
    def __init__(self, target):
        super(PymolListener, self).__init__()
        self.target = target
        self.configure_controller()

    def configure_controller(self):
        self.controller = Leap.Controller()
        self.controller.set_policy(Leap.Controller.POLICY_OPTIMIZE_HMD);
        self.controller.add_listener(self)
        self.configure_gestures()

    def configure_gestures(self):
        self.controller.enable_gesture(Leap.Gesture.TYPE_SWIPE)
        self.controller.enable_gesture(Leap.Gesture.TYPE_CIRCLE)

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

@coroutine
def gesture_dispatch(target_map):
    previous_frame = None
    while True:
        frame = (yield)
        if previous_frame:
            gestures = frame.gestures(previous_frame)
        else:
            gestures = frame.gestures()
        previous_frame = frame
        for gesture in gestures:
            if gesture.is_valid and gesture.type in target_map:
                target = target_map[gesture.type]
                target.send(gesture)

@coroutine
def circle_gesture(target):
    previous_gesture_id = None
    while True:
        gesture = (yield)
        if gesture.id != previous_gesture_id:
            target.send([(0.0, 0.0, 0.0), None])
        circle = Leap.CircleGesture(gesture)
        if (circle.pointable.direction.angle_to(circle.normal) <= Leap.PI/2):
            # clockwise
            target.send([circle.normal, circle.progress])
        else:
            # counterclockwise
            target.send([circle.normal, -circle.progress])
        previous_gesture_id = gesture.id


@coroutine
def scale_z(factor):
    previous_scale = None
    while True:
        vector, scale = (yield)
        print 'scale_z', vector, scale
        if previous_scale and scale and scale != previous_scale:
            delta_scale = (scale - previous_scale) * factor
            print 'delta_scale', delta_scale
            view = list(cmd.get_view())
            view[11] += delta_scale
            view[15] -= delta_scale
            view[16] -= delta_scale
            cmd.set_view(view)
        previous_scale = scale

@coroutine
def rotate(factor):
    previous_rotation = None
    while True:
        vector, rotation = (yield)
        if previous_rotation and rotation and previous_rotation != rotation:
            delta_rotation = (rotation - previous_rotation) * factor
            cmd.rotate(vector.to_float_array(), delta_rotation)

        previous_rotation = rotation

@coroutine
def swipe_gesture(target):
    previous_gesture_id = None
    while True:
        gesture = (yield)
        swipe = Leap.SwipeGesture(gesture)
        distance = swipe.start_position.distance_to(swipe.position)
        print 'swipe', swipe.direction, distance
        target.send([swipe.direction, distance])

cmd.fetch('1rx1')
cmd.orient()

listener = PymolListener(gesture_dispatch({Leap.Gesture.TYPE_SWIPE: swipe_gesture(scale_z(0.5)),
                                           Leap.Gesture.TYPE_CIRCLE: circle_gesture(rotate(45.0))}))
