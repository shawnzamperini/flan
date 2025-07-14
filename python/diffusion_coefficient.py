import numpy as np

# Doesn't totally work as intended yet. This is calculating entire averages over a frame - not taking into account radial differences
def const_x_avg_velocity(frame_number: int, x_coordinate: float):
    imp_vX_array = fp.load_data_frame("imp_vX", frame_number)
    imp_vX_values = imp_vX_array[:, :, 8].flatten()
    imp_vX_avg = [val for val in imp_vX_values if val != 0]
    velAv = np.mean(imp_vX_avg)
    return velAv


