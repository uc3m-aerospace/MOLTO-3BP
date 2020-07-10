def load_kernels():

    import spiceypy as spice
    import os

    dirAppend = "kernels" + os.sep

    spice.furnsh(dirAppend + "naif0010.tls")
    spice.furnsh(dirAppend + "2000001.bsp")
    spice.furnsh(dirAppend + "de432s.bsp")
    spice.furnsh(dirAppend + "gm_de431.tpc")
    spice.furnsh(dirAppend + "pck00010.tpc")
