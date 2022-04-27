from math import sin, cos, sqrt, atan, degrees
import math
import numpy as np
import statistics as st

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
            self.m = 0.999923
            self.m_0 = 0.9993
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2flh(self, X, Y, Z):
        r = sqrt(X**2 + Y**2)
        lat_prev = atan(Z / (r * (1 - self.ecc2)))
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat_prev)) ** 2)
        h = r / cos(lat_prev) - N
        lat_next = atan((Z / r) * (((1 - self.ecc2 * N / (N + h))**(-1))))
        epsilon = 0.0000001 / 206265
        while abs(lat_prev - lat_next) < epsilon:
            lat_prev = lat_next
            N = self.a / sqrt(1 - self.ecc2 * (sin(lat_prev)) ** 2)
            h = r / cos(lat_prev) - N
            lat_next = atan((Z / r) * (((1 - self.ecc2 * N / (N + h))**(-1))))
        phi = lat_prev
        lam = atan(Y / X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(phi)) ** 2)
        h = r / cos(phi) - N
        return degrees(phi), degrees(lam), h
    
    def flh2xyz(self, phi, lam, h):
        N = self.a / sqrt(1 - self.ecc2 * (sin(phi)) ** 2)
        X = (N + h) * cos(phi) * cos(lam)
        Y = (N + h) * cos(phi)*sin(lam)
        Z = (N * (1 - self.ecc2) + h) * sin(phi)
        return X, Y, Z
    
    def u1992(self, phi, lam):
        N = self.a/(np.sqrt(1-self.ecc2 * np.sin(phi)**2))
        t = np.tan(phi)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (np.cos(phi))**2
        
        lam_00 = math.radians(19)
        l = lam - lam_00
        
        A_0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A_2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A_4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A_6 = (35 * (self.ecc2**3))/3072 
        
    
        sigma = self.a * ((A_0*phi) - (A_2*np.sin(2*phi)) + (A_4*np.sin(4*phi)) - (A_6*np.sin(6*phi)))
        
        x = sigma + ((l**2)/2) * N *np.sin(phi) * np.cos(phi) * (1 + ((l**2)/12) * ((math.cos(phi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((math.cos(phi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*math.cos(phi)) * (1 + ((((l**2)/6) * (math.cos(phi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (math.cos(phi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
        x1992 = round(self.m_0*x - 5300000, 3)
        y1992 = round(self.m_0*y + 500000, 3)
        return x1992, y1992 

    def u2000(self, phi, lam):
        N = self.a/math.sqrt(1-self.ecc2*math.sin(phi)**2)
        t = np.tan(phi)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (np.cos(phi))**2
        lam = math.degrees(lam)
        if 16.5 > lam > 13.5 :
            s = 5
            lam0 = 15
        elif 19.5 > lam > 16.5 :
            s = 6
            lam0 = 18
        elif 22.5 > lam > 19.5:
            s = 7
            lam0 = 21
        #elif 25.5 > lam > 22.5:
        else:
            s = 8
            lam0 = 24
        lam = math.radians(lam)
        lam0 = math.radians(lam0)
        l = lam - lam0
    
        A_0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A_2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A_4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A_6 = (35 * (self.ecc2**3))/3072 
        
        sigma = self.a * ((A_0*phi) - (A_2*np.sin(2*phi)) + (A_4*np.sin(4*phi)) - (A_6*np.sin(6*phi)))
        x = sigma + ((l**2)/2) * N *np.sin(phi) * np.cos(phi) * (1 + ((l**2)/12) * ((math.cos(phi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((math.cos(phi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*math.cos(phi)) * (1 + ((((l**2)/6) * (math.cos(phi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (math.cos(phi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
        x2000 = round(self.m * x, 3) 
        y2000 = round(self.m * y + (s*1000000) + 500000, 3)
        return(x2000, y2000)

    def xyz2elewacja(self,X,Y,Z,phi,lam,h):
        Na = self.a/math.sqrt(1 - self.ecc2 *((math.sin(np.deg2rad(phi))**2)))
        XYZs = np.array([X,Y,Z])
    
        Xa = (Na+h)*math.cos(phi)*math.cos(lam)
        Ya = (Na+h)*math.cos(phi)*math.sin(lam)
        Za = (Na*(1-self.ecc2)+h)*math.sin(phi)
        XYZa = np.array([Xa,Ya,Za])
        dX = XYZs-XYZa
        R = np.array([[-sin(phi)*cos(lam), -sin(phi)*sin(lam), cos(phi)],
                      [-sin(lam), cos(lam), 0],
                      [cos(phi)*cos(lam), cos(phi)*sin(lam), sin(phi)]])
    
        neu = R.dot(dX)
    
        Az = math.atan2(neu[1],neu[0])
        if Az<0:
            Az += 2*math.pi
        #print('Az:',np.rad2deg(Az))
        
        el = math.asin(neu[2]/(math.sqrt(neu[0]**2+neu[1]**2+neu[2]**2)))
        #print('elewacja:',np.rad2deg(el))
        return el,Az