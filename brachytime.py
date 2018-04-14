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

import numpy as np
import math as m
from datetime import datetime as dt
import scipy.interpolate as spi

__version__ = '2018-04-13T17:24:00Z'

#Use linear interpolation for wires of intermediate length
def interpdose(wlen, dist):
    lens  = [50,60,70]
    doseatlens = [calcdoserate(i, dist) for i in lens]
    y_interp = spi.interp1d(lens, doseatlens)
    drate = float(y_interp(wlen))
    return drate

#Return the dose rate in Gy/h at a given distance (mm) from a wire of length between 50 mm and 70 mm
#Formulae come from quadratic OLS fits of log10(doserate) to log10(dist)
def calcdoserate(wlen, dist):
    dist = np.log10(dist)
    if wlen == 50:
        drate = -0.3738*dist**2 -0.5991*dist + 0.3876 
    elif wlen == 60:
        drate = -0.3006*dist**2 -0.7100*dist + 0.4438
    elif wlen == 70:
        drate = -0.2891*dist**2 -0.6861*dist + 0.4288 
    elif (50 < wlen < 60) or (60 < wlen < 70):
        drate = np.log10(interpdose(wlen,np.power(10,dist)))
    else:
        raise ValueError('Wire length outside range [50,70]!')
    return np.power(10,drate)

#Return the treatment time in hours
def treattime(sep, ulens, llens, delakr, deltime, midtime, halflife, presc):
    
    #If no halflife is specified, assume it's Ir-192
    if halflife == None:
        halflife = 6379000
    
    #Ensure dates formatted according to ISO 8601
    try:
        deldate=dt.strptime(deltime, '%Y-%m-%dT%H:%M')
        middate=dt.strptime(midtime, '%Y-%m-%dT%H:%M')        
    except Exception as e:
        raise ValueError('Dates and times must be formatted according to ISO 8601 (YYYY-MM-DDTHH:MM)')

    if not (len(ulens) == len(llens)) and not (len(ulens) == len(llens)-1):
        raise ValueError('Numbers of wires in upper row must be equal to or one less than number of wires in lower row.')
    if any(x < 50 or x > 70 for x in ulens+llens):
        raise ValueError('All wire lengths must be between 50 mm and 70 mm.')
    if middate < deldate:
        raise ValueError('Mid-treatment time must be after source delivery time!')
    
    uw_coords = []
    lw_coords = []
    bdp_coords = []
    
    #Compute coordinates of upper and lower wires
    for i in range(len(ulens)):
        uw_coords.append((sep*(i+0.5),(m.sqrt(3)/2*sep)))
    
    for i in range(len(llens)):
        lw_coords.append((sep*i,0))
    
    #Is total number of wires (hence triangles and BDPs) odd or even?
    if (len(llens)+len(ulens))%2 == 1:
        for i in range(len(llens)+1):
            bdp_coords.append((((i+1)/2)*sep,sep/(((i+1)%2+1)*m.sqrt(3))))
            
    elif (len(llens)+len(ulens))%2 == 0:
        for i in range(len(llens)):
            bdp_coords.append((((i+1)/2)*sep,sep/(((i+1)%2+1)*m.sqrt(3))))

    #A list of lists containing the dose rates at each BDP from each wire
    doserates = [[] for _ in range(len(bdp_coords))]

    for ix, i, jx, j in ((ix, i, jx, j) for ix, i in enumerate(bdp_coords) for jx, j in enumerate(uw_coords+lw_coords)):
        doserates[ix].append(calcdoserate((ulens+llens)[jx],m.hypot(i[0]-j[0], i[1]-j[1])))
      
    #Decay-correct the source activity to mid-treatment time, then use this to compute the MBDR and thus treatment time
    midakr = delakr*m.exp(-m.log(2)*((middate-deldate).total_seconds())/halflife)
    mbdr = midakr*np.mean([sum(i) for i in doserates])    
    ttime = presc/(0.85*mbdr) #The RDR is 0.85*MBDR
    
    return float('%.1f' % round(ttime, 1))
 
if __name__ == "__main__":
    print('This is brachytime version {version}'.format(version=__version__))
    exit()
