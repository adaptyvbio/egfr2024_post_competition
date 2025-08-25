from loguru import logger

def set_up_pyrosetta():
    """Initialize PyRosetta with required settings"""
    logger.info("Setting up PyRosetta")
    import pyrosettacolabsetup
    pyrosettacolabsetup.install_pyrosetta(
        serialization=True, 
        cache_wheel_on_google_drive=False
    )
    logger.success("PyRosetta setup complete")

def initialize_pyrosetta():
    """Initialize PyRosetta with specific options"""
    import pyrosetta as pr
    logger.info("Initializing PyRosetta")
    pr.init(extra_options=f"-ignore_unrecognized_res "
            "-ignore_zero_occupancy "
            "-mute all "
            "-holes:dalphaball /root/bindcraft/functions/DAlphaBall.gcc "
            "-corrections::beta_nov16 true "
            "-relax:default_repeats 1 "
            )
    logger.info("PyRosetta initialized successfully") 