# ecef_to_eci.py

# Usage: python3 ecef_to_eci.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km
# Description: Converts Year, Month, Day, Hour, Minute, Second, ecef_x_km, ecef_y_km, ecef_z_km inputs from the ECEF frame to ECI coordinates.

# Parameters:
#  year
#  month
#  day
#  hour
#  minute
#  second
#  ecef_x_km
#  ecef_y_km
#  ecef_z_km

# Output:
#  Coordinates in the eci reference frame (x, y, z), answer in km

# Written by Owen Davies
# Other contributors: None

# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
import math
import sys

# "constants"
R_E_KM = 6378.137
e_E = 0.081819221456
r_E_km = 6378.1363
w = 0.00007292115
Delta_UT1 = 236.555

# helper functions - Denominator calculation
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0-ecc*ecc*math.pow(math.sin(lat_rad),2.0))

# initialize script arguments
year: float(sys.argv[1])
month: float(sys.argv[2])
day: float(sys.argv[3])
hour: float(sys.argv[4])
minute: float(sys.argv[5])
second: float(sys.argv[6])
ecef_x_km: float(sys.argv[7])
ecef_y_km: float(sys.argv[8])
ecef_z_km: float(sys.argv[9])

# parse script arguments
if len(sys.argv)==10:
  year = float(sys.argv[1])
  month = float(sys.argv[2])
  day = float(sys.argv[3])
  hour = float(sys.argv[4])
  minute = float(sys.argv[5])
  second = float(sys.argv[6])
  ecef_x_km = float(sys.argv[7])
  ecef_y_km = float(sys.argv[8])
  ecef_z_km = float(sys.argv[9])
else:
  print(\
   'Usage: '\
   'python3 eci_ecef.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km'\
  )
  exit()
  
# Calculate JD_frac
jd = day-32075+(1461*(year+4800+(month-14)/12)/4)+(367*(month-2-((month-14)/12)*12)/12)-(3*((year+4900+((month-14)/12))/100)/4)
jdint = int(jd)
jdmidnight = jdint-0.5
Dfrac = (second+60*(minute+60*hour))/86400
JDfrac = jdmidnight+Dfrac 

# Calculate GMST Angle
Tut1 = (JDfrac - 2451545.0)/36525
GMSTthetaseconds = 67310.54841+(876600*60*60+8640184.812866)*Tut1+0.093104*(Tut1**2)+(-6.2*(10**-6))*(Tut1**3)
GMSTangrad = (GMSTthetaseconds % 86400)*w + 2*math.pi
GMSTrad = (math.fmod(GMSTangrad, 2*math.pi))

# Calculate the eci vector, use the transpose of the Rz matrix * ECEF vector
eci_x_km = (ecef_x_km*math.cos(-GMSTrad))+(ecef_y_km*(math.sin(-GMSTrad)))
eci_y_km = (ecef_x_km*-math.sin(-GMSTrad))+(ecef_y_km*math.cos(-GMSTrad))
eci_z_km = (ecef_z_km)

# Display final answers
print(f"{eci_x_km:.3f}")
print(f"{eci_y_km:.3f}")
print(f"{eci_z_km:.3f}")

