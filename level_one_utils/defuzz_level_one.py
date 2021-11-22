"""
This class defuzzifies the Overlap and Spatial Relationship fuzzy memberships given crisp input values
"""

import numpy as np
import skfuzzy as fuzz
from skfuzzy import control as ctrl


class Defuzz:
    def __init__(self):
        self.overlap = ctrl.Antecedent(universe=np.arange(-1.1, 1.1, 0.1), label='overlap')
        self.sr = ctrl.Antecedent(universe=np.arange(-1, 361, 1), label='spatial_relationships')

        self.overlap['Overlap'] = fuzz.trapmf(self.overlap.universe, [0.0, 0.2, 0.7, 1.0])
        self.overlap['No Overlap'] = fuzz.trapmf(self.overlap.universe, [-1.0, -0.7, -0.2, 0.0])

        # 0 < HOF < 30 | 331 < HOF < 360: Right
        self.sr['Right1'] = fuzz.trimf(self.sr.universe, [-1, 15, 31])
        self.sr['Right2'] = fuzz.trimf(self.sr.universe, [330, 345, 360])
        # 31 < HOF < 60: Above Right
        self.sr['Above Right'] = fuzz.trimf(self.sr.universe, [30, 45, 61])
        # 61 < HOF < 120: Above
        self.sr['Above'] = fuzz.trimf(self.sr.universe, [60, 90, 121])
        # 121 < HOF < 150: Above Left
        self.sr['Above Left'] = fuzz.trimf(self.sr.universe, [120, 135, 151])
        # 151 < HOF < 210: Left
        self.sr['Left'] = fuzz.trimf(self.sr.universe, [150, 180, 211])
        # 211 < HOF < 240: Below Left
        self.sr['Below Left'] = fuzz.trimf(self.sr.universe, [210, 225, 241])
        # 241 < HOF < 300: Below
        self.sr['Below'] = fuzz.trimf(self.sr.universe, [240, 270, 301])
        # 301 < HOF < 330: Below Right
        self.sr['Below Right'] = fuzz.trimf(self.sr.universe, [300, 315, 331])

    def defuzzify_results(self, iou, sr_angle):
        # Compute the overlap result
        overlap = fuzz.interp_membership(self.overlap.universe, self.overlap['Overlap'].mf, iou)
        no_overlap = fuzz.interp_membership(self.overlap.universe, self.overlap['No Overlap'].mf, iou)
        membership = {'Overlap': overlap, 'No Overlap': no_overlap}
        overlap_label = max(membership, key=membership.get)

        # Compute the spatial relationship result
        right_1 = fuzz.interp_membership(self.sr.universe, self.sr['Right1'].mf, sr_angle)
        right_2 = fuzz.interp_membership(self.sr.universe, self.sr['Right2'].mf, sr_angle)
        above_right = fuzz.interp_membership(self.sr.universe, self.sr['Above Right'].mf, sr_angle)
        above = fuzz.interp_membership(self.sr.universe, self.sr['Above'].mf, sr_angle)
        above_left = fuzz.interp_membership(self.sr.universe, self.sr['Above Left'].mf, sr_angle)
        left = fuzz.interp_membership(self.sr.universe, self.sr['Left'].mf, sr_angle)
        below_left = fuzz.interp_membership(self.sr.universe, self.sr['Below Left'].mf, sr_angle)
        below = fuzz.interp_membership(self.sr.universe, self.sr['Below'].mf, sr_angle)
        below_right = fuzz.interp_membership(self.sr.universe, self.sr['Below Right'].mf, sr_angle)
        membership = {'Right1': right_1, 'Right2': right_2, 'Above Right': above_right, 'Above': above,
                      'Above Left': above_left, 'Left': left, 'Below Left': below_left, 'Below': below,
                      'Below Right': below_right}
        sr_label = max(membership, key=membership.get)
        if sr_label == 'Right1' or sr_label == 'Right2':
            sr_label = 'Right'
        return overlap_label, sr_label
