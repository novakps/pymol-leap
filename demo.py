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

class PymolListener(Leap.Listener):
    def __init__(self, target):
        super(PymolListener, self).__init__()
        self.target = target
        self.configure_controller()

    def configure_controller(self):
        self.controller = Leap.Controller()
        self.controller.set_policy(Leap.Controller.POLICY_OPTIMIZE_HMD);
        self.controller.add_listener(self)
        # self.configure_gestures()

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
def frame_scaler():
    previous_frame = None
    while True:
        frame = (yield)
        if previous_frame and frame and frame.scale_probability(previous_frame) > 0.2:
            scale_factor = frame.scale_factor(previous_frame)
            if scale_factor != 1.0:
                view = list(cmd.get_view())
                delta_z = math.log(scale_factor) * 100.0
                view[11] += delta_z
                view[15] -= delta_z
                view[16] -= delta_z
                cmd.set_view(view)
        previous_frame = frame

@coroutine
def hand_rotator():
    previous_frame = None
    while True:
        hand = (yield)
        if previous_frame and hand and hand.rotation_probability(previous_frame) > 0.1:
            view = list(cmd.get_view())
            matrix = hand.rotation_matrix(previous_frame)
            matrix *= Leap.Matrix(Leap.Vector(*view[0:3]),
                                  Leap.Vector(*view[3:6]),
                                  Leap.Vector(*view[6:9]))
            view[:9] = matrix.to_array_3x3()
            cmd.set_view(view)
        if hand:
            previous_frame = hand.frame
        else:
            previous_frame = None

@coroutine
def hands_splitter(target):
    while True:
        frame = (yield)
        for hand in frame.hands:
            target.send(hand)

@coroutine
def broadcaster(targets):
    while True:
        frame = (yield)
        for target in targets:
            target.send(frame)

@coroutine
def multi_hands_filter(target):
    while True:
        frame = (yield)
        if frame and len(list(frame.hands)) > 1:
            target.send(frame)
        else:
            target.send(None)

@coroutine
def pinch_filter(target):
    while True:
        hand = (yield)
        if hand.pinch_strength > 0.5:
            target.send(hand)
        else:
            target.send(None)

@coroutine
def select_all():
    previous_hand = None
    while True:
        hand = (yield)
        if hand and not previous_hand:
            cmd.select('all')
        elif previous_hand and not hand:
            cmd.select('none')
        previous_hand = hand

cmd.fetch('1rx1')
cmd.orient()

listener = PymolListener(
    broadcaster([
        hands_splitter(
            pinch_filter(
                broadcaster([
                    hand_rotator(),
                    select_all()]
                ))
        ),
        multi_hands_filter(
            frame_scaler()
        )
    ]))
