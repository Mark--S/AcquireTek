import argparse
import scope_connections
import scopes
import utils
from datetime import datetime

def save_scopeTraces(fileName, scope, channel):
    """Save a number of scope traces to file - uses compressed .pkl"""

    datestr = datetime.now().strftime("%d-%b_%H:%M:%S")
    scope._get_preamble(channel)
    results = utils.PickleFile("%s_%s.pkl" % (fileName, datestr), 1)
    results.add_meta_data("timeform_1", scope.get_timeform(channel))

    t_start, loopStart = time.time(),time.time()
    while(time.time()-t_start < 300): # Run for 5 mins
        try:
            results.add_data(scope.get_waveform(channel), 1)
        except Exception, e:
            print "Scope died, acquisition lost."
            print e
            results.save()
            results.close()
        if i % 100 == 0 and i > 0:
            print "%d traces collected" % (i)
            loopStart = time.time()
    print "%d traces collected TOTAL - took : %1.1f s"a % (i, (time.time()-t_start))
    results.save()
    results.close()
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="channel", type=int, 
                        default=1,
                        help="Scope channel to acquire / trigger from")
    parser.add_argument("-t", dest="trigger", type=float,
                        help="Trigger level (if positive / negative will trigger off rising / falling edge")
    parser.add_argument("-n", dest="name", type=str,
                        default="test",
                        help = "Base name for .pkl file that traces will be saved too")
    args = parser.parse_args()    
    
    ##########################################
    usb_conn = scope_connections.VisaUSB()
    scope = scopes.Tektronix3000(usb_conn)    
    ##########################################
    scope_chan = args.channel # Set channel to be the one passed
    if args.trigger < 0:
        falling_edge = True
    else:
        falling_edge = False
    trigger = args.trigger
    termination = 50 # Ohms
    y_div_units = 0.2 # volts
    x_div_units = 10e-9 # seconds
    y_offset = -2.5*y_div_units # offset in y (2.5 divisions up)
    #y_offset = 0.5*y_div_units # offset in y (for UK scope)
    x_offset = +2*x_div_units # offset in x (2 divisions to the left)
    record_length = 1e3 # trace is 1e3 samples long
    half_length = record_length / 2. # For selecting region about trigger point
    ###########################################
    scope.unlock()
    scope.set_horizontal_scale(x_div_units)
    scope.set_horizontal_delay(x_offset) #shift to the left 2 units
    scope.set_channel_y(scope_chan, y_div_units, pos=2.5)
    scope.set_channel_termination(scope_chan, termination)
    scope.set_single_acquisition() # Single signal acquisition mode
    scope.set_record_length(record_length)
    scope.set_data_mode(half_length-80, half_length+20)
    #scope.set_edge_trigger(-0.5*y_div_units, scope_chan, falling=True)
    scope.set_edge_trigger(trigger, scope_chan, falling=falling_edge)
    scope.lock()
    scope.begin() # Acquires the pre-amble! 
