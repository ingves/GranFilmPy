
import numpy as np


# --------------- Function to run the GranFilm executable --------------- #


def run_command( cmd, wd ):
    '''
    Runs the given command in the given directory.
    Returns:
    stdout
    stderr
    returncode
    time
    '''
    import os
    import time
    from subprocess import Popen, PIPE
    #
    curr_path = os.getcwd()
    os.chdir(wd)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    start = time.time()
    (out, err) = p.communicate()
    runningtime = time.time() - start
    os.chdir(curr_path)
    return (out, err, p.returncode, runningtime)


# --------------- Get data from file and return them as an array--------------- #


def get_results_from_file( filename ):
    '''
    Return the content of the file filename and returns 
    the content as numpy array
    '''
    try:
        data = np.genfromtxt(filename, comments='#')
    except:
        print ("ERROR : reading data file %s" % filename)
    #    
    return data

      
# --------------- Time info --------------- #


def timestamp():
    '''
    Returns a time stamp string
    '''
    try:
        # make the time stamp string
        import time
        import datetime  
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    except:
        print ("ERROR : making the time-stamp-string")
    #    
    return st


# --------------- Data format --------------- #


def reformat(data,fmt,deg="True"):
    '''
    Reformates the output from an numpy.array

    Options
    -------
    data : numpy.array
       Data to be reformated
    fmt  : string
       Format to use for the output
       Possible values 'Complex', 'Real', 'Imag', 'Amplitude', and 'Phase'
    deg  :  string 
       If phase is returned, use degreee or radians
       Default : Degree
    '''
    import numpy
    # Check type of data
    FMT={'Complex', 'Real', 'Imag', 'Amplitude', 'Phase' }
    assert (fmt in FMT), "ERROR : Unsupported format %s " % fmt


    if (fmt=='Complex'):
        return data
    elif (fmt=='Real'): 
        return data.real
    elif (fmt=='Imag'):
        return data.imag
    elif (fmt=='Amplitude'):
        return np.abs(data)
    elif (fmt=='Phase'):
        return np.angle( data, deg=deg )
