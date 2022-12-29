import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import readEqdsk
import KS2D
from matplotlib import animation
from scipy.stats.kde import gaussian_kde
import seaborn as sns
from shapely.geometry import Polygon, Point
import geopandas as gpd
from scipy.spatial import cKDTree as KDTree

from out_functions import *
from glob_vars import *

# =============================================================================
#     normalize units
# =============================================================================
    checkvals(dist0)
    checkvals(dist1)
    checkvals(traj)
    # checkvals(poink)
    
# =============================================================================
#     dist0 plotting
# =============================================================================
    # dist0.hist(bins = 50, figsize= (9,6))
    # tok_scatter(dist0,'xflr','zflr','gist_rainbow')
    # tok_scatter(dist0,'xgc','zgc','Set1')
    # tok_scatter_otp(dist0,'xgc','zgc')
    
# =============================================================================
#     dist1 plotting
# =============================================================================
    # dist1.hist(bins = 50, figsize= (9,6))
    # tok_scatter(dist1,'xgc','zgc','Set1')
    # tok_scatter(dist1,'xflr','zflr','Set1')
    # tok_scatter_otp(dist1,'xgc','zgc')
    # tok_scatter_otp(dist1,'xflr','zflr')
    # heatmap(dist1, 'xgc', 'zgc',1)
    # heatmap(dist1, 'xflr', 'zflr',1)
    # heatmap(dist1, 'xgc', 'zgc',0)
    # snsheatmap(dist1, 'xflr', 'zflr')
    enptchheatmap(dist1, 'ptch', 'en', 1)
    
# =============================================================================
#     trajectory plotting
# =============================================================================
    # trajTok(traj,False)
    # trajside(traj)
    
# =============================================================================
#     poincare plotting
# =============================================================================
    # poink.plot.scatter('xgc','pol')
    # poink.plot.scatter('pol','thet')
    
# =============================================================================
#     other functions
# =============================================================================
    # points_out(dist1, 'xflr', 'zflr')
    # KS_test(dist1)
    # KLdivergence(dist1)
    # gc_flrdiff(traj)
    # plt.hist2d(dist0.rhol,dist0.en,bins=50)