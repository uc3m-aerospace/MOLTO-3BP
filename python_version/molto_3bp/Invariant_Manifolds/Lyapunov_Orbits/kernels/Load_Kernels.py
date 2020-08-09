def load_kernels():

    import spiceypy as spice
    import os

    dirAppend = "Lyapunov_Orbits" + os.sep + "kernels" + os.sep

    spice.furnsh(dirAppend + "gm_de431.tpc")
