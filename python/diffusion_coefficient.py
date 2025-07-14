import numpy as np
import flan_plots
fp=flan_plots.FlanPlots(self.nc)

# Takes a frame, radial index, and physical quantity. Calculates the average of the physical quantity in a given frame at that radial index. 
def phAvg(int radialIndex, int frame, string physicalQuantity) -> float:
    a, b, c = fp.plot_frame_xy(physicalQuantity, frame=frame, z0=0.01, norm_type="log", min=1e-5, max=1e-1)
    
    total = 0.0
    for i in range(96):  # i from 0 to 95
        if(c[i][radialIndex != 0) 
            total += c[i][radialIndex]
    
    phAvg = total / 96.0
    return phAvg






