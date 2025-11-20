import Observation
import glob, os
package_path = os.path.dirname(os.path.dirname(os.path.abspath(Observation.__file__)))
from Observation import *




data_path = os.path.join(package_path, "data")
data_cube__path = os.path.join(data_path, "Emission_cube")
if os.path.exists(data_cube__path) is False:
    os.makedirs(data_cube__path)
for name in [ "cube_01.fits", "lya_cube_merged_with_artificial_source_CU_1pc.fits", "CGM_cube.fits", "galaxy_disk_cube.fits", "galaxy_and_cgm_cube.fits"]:
    if os.path.exists(os.path.join(data_cube__path, "cube_01.fits")) is False:
        a = download(
            url="https://nuage.osupytheas.fr/s/fo4AKjoTZ4fBytw/download/"+name,
            file=os.path.join(data_cube__path, name))

instruments, database = load_instruments(sheet_id="1Ox0uxEm2TfgzYA6ivkTpU4xrmN5vO5kmnUPdCSt73uU",sheet_name="instruments.csv",database="Online DB")


FB = ExposureTimeCalulator(instrument="SCWI PERF",instruments=instruments,database=database, Δλ=-10,Signal=5e-17,lambda_stack=10,interpolation="gaussian",Line_width=6)#, spectra="↳ cube_01-resampled_phys",min=0.54,max=0.85
# FB = ExposureTimeCalulator(instrument="SCWI MCP",instruments=instruments,database=database, Δλ=-10,Signal=5e-18,lambda_stack=10, spectra="↳ cube_01-resampled_phys",min=0.54,max=0.85,interpolation="gaussian",sky_lines=False,atmlambda=False,QElambda=False)
# FB = ExposureTimeCalulator(instrument="SCWI PERF",instruments=instruments,database=database)
