from mcmax3d_analysis.mcmax3d_observables import convert_flux
from mcmax3d_analysis.mcmax3d_image import display_image
from mcmax3d_analysis.mcmax3d_profiles import radial_brightness_profile
import numpy as np
import matplotlib.pyplot as plt
import sys
plt.style.use('fancy')

print(radial_brightness_profile('/data/users/bportilla/runs/MCMax3D_images/run_88','I'))

