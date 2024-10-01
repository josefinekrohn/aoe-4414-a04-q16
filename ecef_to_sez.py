# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  Converts ECEF vector components to SEZ
# Parameters:
#  o_x_km: origin ECEF x-component in km
#  o_y_km: origin ECEF y-component in km
#  o_z_km: origin ECEF z-component in km
#  x_km: ECEF x-component in km
#  y_km: ECEF x-component in km
#  z_km: ECEF x-component in km
#  ...
# Output:
#  #  Prints the SEZ vector (south vector, east vector, and vector normal to the ellipsoid surface) in km
#  Inputs an origin ECEF location and a position vector in ECEF and prints out position in SEZ-coordinates
#
# Written by Josefine Krohn
# Other contributors: Brad Denby (format & ecef_to_llh.py) 
#

# import Python modules
import math # math module
import sys # argv

# "constants"
R_E_KM = 6378.1363
E_E = 0.081819221456  

# helper functions

## calculated denominator
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

# Matrix Multiplication
def matrix_multiplication(A, B):
    result = [[0 for column_result in range(len(B[0]))] for row_result in range(len(A))] # setting up results vector
    for row_A in range(len(A)): # loop for rows of matrix A
        a = A[row_A] # creating a vector of the current row of A
        for column_B in range(len(B[0])): # loop for columns of matrix B
            b = [row[column_B] for row in B] # creating a vector of the current column of B       
            result_value = 0 # initializing result value
            for column_A in range(len(b)): # loop for # of values in vector b (equals # of values in a)
                value = a[column_A] * b[column_A] # multiplying values in corresponding rows and columns
                result_value = result_value + value # adds the row*column values together
        result[row_A][column_B] = result_value # updates result matrix
    return result  

# Converts ECEF vector components to LLH 
def ecef_to_llh(r_x_km, r_y_km, r_z_km):  
    # calculate longitude
    lon_rad = math.atan2(r_y_km,r_x_km)
    # initialize lat_rad, r_lon_km, r_z_km
    lat_rad = math.asin(r_z_km/math.sqrt(r_x_km**2+r_y_km**2+r_z_km**2))
    r_lon_km = math.sqrt(r_x_km**2+r_y_km**2)
    prev_lat_rad = float('nan')
    # iteratively find latitude
    c_E = float('nan')
    count = 0
    while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
        denom = calc_denom(E_E,lat_rad)
        c_E = R_E_KM/denom
        prev_lat_rad = lat_rad
        lat_rad = math.atan((r_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
        count = count+1
    # calculate hae
    hae_km = r_lon_km/math.cos(lat_rad)-c_E
    return lat_rad, lon_rad, hae_km # latitude (rad), longitude (rad), and HAE (km)



# initialize script arguments
o_x_km = float('nan') # origin ECEF x-component in km
o_y_km = float('nan') # origin ECEF y-component in km
o_z_km = float('nan') # origin ECEF z-component in km
x_km = float('nan') # ECEF x-component in km
y_km = float('nan') # ECEF y-component in km
z_km = float('nan') # ECEF z-component in km

# parse script arguments
if len(sys.argv)==7:
  o_x_km = float(sys.argv[1])
  o_y_km = float(sys.argv[2])
  o_z_km = float(sys.argv[3])
  x_km = float(sys.argv[4])
  y_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
   'Usage: '\
   'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
  )
  exit()





# write script below this line

# determine ECEF vector from the station to the object
r_x = x_km - o_x_km
r_y = y_km - o_y_km
r_z = z_km - o_z_km

r_sez = [[r_x], [r_y], [r_z]] # SEZ vector pre-rotation


# converting station ECEF vector to LLH
o_lat_rad, o_lon_rad, o_hae_km = ecef_to_llh(o_x_km, o_y_km, o_z_km)


# SEZ to ECEF rotation
phi = o_lat_rad
theta = o_lon_rad
R_y_inv = [[math.cos(theta), math.sin(theta), 0], [-math.sin(theta), math.cos(theta), 0], [0, 0, 1]]
R_z_inv = [[math.sin(phi), 0, -math.cos(phi)], [0, 1, 0], [math.cos(phi), 0, math.sin(phi)]]
R_y_inv_r_sez = matrix_multiplication(R_y_inv,r_sez) # R_y(90deg - phi)*r_sez
r_sez = matrix_multiplication(R_z_inv,R_y_inv_r_sez) # r_ECEF = R_z(theta)*R_y(90deg - phi)*r_sez


s_km = r_sez[0][0]
e_km = r_sez[1][0]
z_km = r_sez[2][0]


print(s_km)
print(e_km)
print(z_km)