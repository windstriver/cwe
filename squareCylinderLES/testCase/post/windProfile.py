# Mean velocity
Uavref = 10    # m/s
zref = 0.364    # m
alpha = 0.326


def Uav(z):
    return Uavref * (z/zref)**(alpha)


# Turbulent intensity
Iuref = 0.208
Ivref = 0.182
Iwref = 0.152
du = 0.191
dv = 0.123
dw = 0.005


def Iu(z):
    return Iuref * (z/zref)**(-du)


def Iv(z):
    return Ivref * (z/zref)**(-dv)


def Iw(z):
    return Iwref * (z/zref)**(-dw)


# Turbulence length scale
Luref = 0.302    # m
Lvref = 0.0815
Lwref = 0.0326
epsu = 0.473
epsv = 0.881
epsw = 1.539
zrefL = 0.254


def Lu(z):
    return Luref * (z/zrefL)**(epsu)


def Lv(z):
    return Lvref * (z/zrefL)**(epsv)


def Lw(z):
    return Lwref * (z/zrefL)**(epsw)


# von Karmon turbulent spectra
def Su(z, f):
    return 4*(Iu(z)*Uav(z))**2*(Lu(z)/Uav(z)) / \
        (1+70.8*(f*Lu(z)/Uav(z))**2)**(5/6)


def Sv(z, f):
    return 4*(Iv(z)*Uav(z)**2*(Lv(z)/Uav(z)) *
              (1+188.4*(2*f*Lv(z)/Uav(z)))**2) / \
              (1+70.8*(f*Lv(z)/Uav(z))**2)**(11/6)


def Sw(z, f):
    return 4*(Iw(z)*Uav(z)**2*(Lw(z)/Uav(z)) *
              (1+188.4*(2*f*Lw(z)/Uav(z)))**2) / \
              (1+70.8*(f*Lw(z)/Uav(z))**2)**(11/6)
