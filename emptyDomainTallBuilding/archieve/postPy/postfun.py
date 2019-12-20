import numpy as np

def integralTimeScale(x, deltaT):
    '''
    Estimate the integral time scale by autocorrelation function
    Input:
        x         signal
        deltaT    time step of the signal
    Output:
        integral time scale
    '''
    N = len(x)
    x = x - np.mean(x)
    autoCorrFun = np.correlate(x, x, mode='full') / np.correlate(x, x)
    mask = np.ones(len(autoCorrFun), dtype=bool)
    mask[:N-1] = False
    autoCorrFun = autoCorrFun[mask]
    M = np.argmax(autoCorrFun < 0) 
    if M == 0:
        return np.trapz(autoCorrFun) * deltaT
    else:
        return np.trapz(autoCorrFun[:M]) * deltaT

