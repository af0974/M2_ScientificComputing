import numpy as np

def def_timewindow(t_start, t_end, t_samp):
    """
    this function returns a time window starting at t_start, finishing at t_end
    with one point every t_samp year
    """
    mynum = int( (t_end - t_start) / t_samp ) 
    twindow = np.linspace(t_start, t_end, num=mynum, endpoint = True, dtype=float)
    return twindow

def cos_tseries(times, period, amp, phase): 
    """
    this function returns amp * cos( 2pi*times/period + phase )
    """
    tseries = amp * np.cos( 2.*np.pi*times/period + phase )
    return tseries 
