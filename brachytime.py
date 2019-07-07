# -*- coding: utf-8 -*-
"""
################################################################################
#                                                                              #
# brachytime                                                                   #
#                                                                              #
################################################################################
#                                                                              #
# LICENCE INFORMATION                                                          #
#                                                                              #
# This program provides a simple means of calculating the treatment time for   #
# an interstitial brachytherapy implant with a geometry comprising equilateral #
# triangles.                                                                   #
#                                                                              #
# (C) 2018 Gavin Donald Kirby                                                  #
# Intial version created 2018-04-11T12:39:46Z                                  #
#                                                                              #
# This software is released under the terms of the GNU General Public License  #
# version 3 (GPLv3).                                                           #
#                                                                              #
# This program is free software: you can redistribute it and/or modify it      #
# under the terms of the GNU General Public License as published by the Free   #
# Software Foundation, either version 3 of the License, or (at your option)    #
# any later version.                                                           #
#                                                                              #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for     #
# more details.                                                                #
#                                                                              #
# For a copy of the GNU General Public License, see                            #
# <http://www.gnu.org/licenses/>.                                              #
#                                                                              #
################################################################################

"""

from datetime import datetime as dt
import math as m
from statistics import mean

import scipy.interpolate as interp

__version__ = '2019-07-07T22:00:00Z'


def calculate_treatment_time(
    separation, upper_lengths, lower_lengths, delivery_akr,
    delivery_datetime, mid_treatment_datetime,
    prescribed_dose, halflife=6379000, printout=False):
    """
    Calculate the treatment time in hours for an applicator of given geometry.
    Applicator separation and upper and lower wire lengths are in millimetres.
    Datetimes are formatted according to ISO-8601, with the exception that
    the seconds and the final "Z" are suppressed, e.g. "2019-07-07T22:00".
    The prescribed dose is in grays, and the AKR is in milligrays per hour.
    If no halflife is specified, default to Ir-192.
    """    
    # Parse datetimes, assuming they are formatted according to ISO 8601.
    try:
        delivery_datetime = dt.strptime(delivery_datetime, '%Y-%m-%dT%H:%M')
        mid_treatment_datetime = dt.strptime(mid_treatment_datetime, '%Y-%m-%dT%H:%M')        
    except Exception:
        raise ValueError('Dates and times must be formatted as YYYY-MM-DDTHH:MM!')

    # Without loss of generality, define the lower row to be the larger.
    if not (len(upper_lengths) == len(lower_lengths)) and not (len(upper_lengths) == len(lower_lengths)-1):
        raise ValueError('Numbers of wires in upper row must be equal to or one less than number of wires in lower row.')
    if any(x < 50 or x > 70 for x in upper_lengths + lower_lengths):
        raise ValueError('All wire lengths must be between 50 mm and 70 mm.')
    if mid_treatment_datetime < delivery_datetime:
        raise ValueError('Mid-treatment time must be after source delivery time!')
    
    # Initialise empty lists to contain coordinates
    # of upper and lower wires, as well as basal dose points.
    upper_coords, lower_coords, bdp_coords = [], [], []
    
    # Compute coordinates of upper and lower wires.
    for i in range(len(upper_lengths)):
        upper_coords.append((separation*(i+0.5), (m.sqrt(3) / 2*separation)))
    
    for i in range(len(lower_lengths)):
        lower_coords.append((separation*i, 0))
    
    # Is total number of wires (hence triangles and BDPs) odd or even?
    if (len(lower_lengths) + len(upper_lengths))%2 == 1:
        for i in range(len(lower_lengths) + 1):
            bdp_coords.append((((i+1)/2)*separation, separation/(((i+1)%2+1)*m.sqrt(3))))
            
    elif (len(lower_lengths)+len(upper_lengths))%2 == 0:
        for i in range(len(lower_lengths)):
            bdp_coords.append((((i+1)/2)*separation, separation/(((i+1)%2+1)*m.sqrt(3))))

    # Populate a list of lists containing the dose rates at each BDP from each wire.
    doserates = [[] for _ in range(len(bdp_coords))]
    for ix, i, jx, j in ((ix, i, jx, j) for ix, i in enumerate(bdp_coords) for jx, j in enumerate(upper_coords + lower_coords)):
        doserates[ix].append(calculate_dose_rate((upper_lengths + lower_lengths)[jx], m.hypot(i[0]-j[0], i[1]-j[1])))
      
    # Decay-correct the source activity to mid-treatment time,
    # then use this to compute the mean basal dose rate and thus treatment time.
    mid_treatment_akr = delivery_akr * m.exp(-m.log(2)*((mid_treatment_datetime - delivery_datetime).total_seconds()) / halflife)
    mbdr = mid_treatment_akr * mean([sum(i) for i in doserates])
    treatment_time = prescribed_dose / (0.85*mbdr)  # The RDR is 0.85*MBDR.

    if not printout:
        return treatment_time
    else:
        print(f'The treatment time is {round(treatment_time, 1)} hours.')


def calculate_dose_rate(wire_length, distance):
    """
    Return the dose rate in Gy/h at a given distance (mm) from a wire of length between 50 mm and 70 mm.
    Formulae come from quadratic OLS fits of log10(doserate) to log10(distance).
    """
    distance = m.log10(distance)
    # Case analysis of wire lengths    
    if wire_length == 50:
        doserate = -0.3738*distance**2 - 0.5991*distance + 0.3876 
    elif wire_length == 60:
        doserate = -0.3006*distance**2 - 0.7100*distance + 0.4438
    elif wire_length == 70:
        doserate = -0.2891*distance**2 - 0.6861*distance + 0.4288 
    elif (50 < wire_length < 60) or (60 < wire_length < 70):
        doserate = m.log10(interpolate_doserate(wire_length, pow(10, distance)))
    else:
        raise ValueError('Wire length outside range [50, 70]!')

    return pow(10, doserate)


def interpolate_doserate(wire_length, distance):
    """Use linear interpolation to calculate dose rates for wires of intermediate length."""
    reference_lengths  = [50, 60, 70]
    doses_at_lengths = [calculate_dose_rate(i, distance) for i in reference_lengths]
    interpolator = interp.interp1d(reference_lengths, doses_at_lengths)

    return float(interpolator(wire_length))  # This is the interpolated doserate.  


if __name__ == "__main__":
    print(f'This is brachytime version {__version__}')
    exit()
