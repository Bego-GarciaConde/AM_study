# -*- coding: utf-8 -*-
#
"""
Created on 29/04/2022
@author: B. Garcia-Conde
"""
import gc
import yt
from yt import YTArray


import numpy as np
from yt.units import G
import array
import pandas as pd

import gc


from config import *
from angular_momentum_comparison import * 
from snapshot_definition import *


def main() :
    snapshot_pre=Snapshot(520)
    for name in snapshots_analysis[1:]:
        snapshot = Snapshot(name)
        angular_momentum_comparison (snapshot, snapshot_pre)
        snapshot_pre = snapshot


if __name__ == "__main__":
    main()