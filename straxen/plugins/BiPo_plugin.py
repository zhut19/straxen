st = straxen.contexts.xenon1t_dali()
@export
class RadonVariables(strax.LoopPlugin):
    __version__ = '0.1.18'
    depends_on = ('events',
                  'event_basics',
                  'peak_basics',
                  'peak_positions',
                  'peak_proximity')
    
    provides = 'Radon_variables'
    def infer_dtype(self):
        dtype = [(('Number of peaks in the event',
                   'n_peaks'), np.int32),]         
        for i in [1, 2]:
            dtype += [
                #main Signal
                ((f'Main S{i} peak index',
                  f's{i}_index'), np.int32),
                ((f'Main S{i} time since unix epoch [ns]',
                  f's{i}_time_0'), np.int64),
                ((f'Main S{i} weighted center time since unix epoch [ns]',
                  f's{i}_center_time_0'), np.int64),
                ((f'Main S{i} area, uncorrected [PE]',
                  f's{i}_area_0'), np.float32),
                #second biggest
                ((f'Alternate S{i} time since unix epoch [ns]',
                  f's{i}_time_1'), np.int64),
                ((f'Alternate S{i} weighted center time since unix epoch [ns]',
                  f's{i}_center_time_1'), np.int64),
                ((f'Area of alternate S{i} in event [PE]',
                  f's{i}_area_1'), np.float32),
                #third biggest
                ((f'third S{i} time since unix epoch [ns]',
                  f's{i}_time_2'), np.int64),
                ((f'third S{i} weighted center time since unix epoch [ns]',
                  f's{i}_center_time_2'), np.int64),
                ((f'Area of third alternate S{i} in event [PE]',
                  f's{i}_area_2'), np.float32),
                ]         
        dtype += [ 
                (('Main S2 reconstructed X position, uncorrected [cm]',
                  's2_x_0'), np.float32),
                (('Main S2 reconstructed Y position, uncorrected [cm]',
                  's2_y_0'), np.float32),
               (('other S2 reconstructed X position, uncorrected [cm]',
                  's2_x_1'), np.float32),
                (('other S2 reconstructed Y position, uncorrected [cm]',
                  's2_y_1'), np.float32),        
                (('second other S2 reconstructed X position, uncorrected [cm]',
                  's2_x_2'), np.float32),
                (('second other S2 reconstructed Y position, uncorrected [cm]',
                  's2_y_2'), np.float32),
                  #fourth biggest
                (('fourth S2 time since unix epoch [ns]',
                  's2_time_3'), np.int64),
                (('fourth S2 weighted center time since unix epoch [ns]',
                  's2_center_time_3'), np.int64),
                (('Area of fourth alternate S2 in event [PE]',
                  's2_area_3'), np.float32),
                (('third other S2 reconstructed X position, uncorrected [cm]',
                  's2_x_3'), np.float32),
                (('third other S2 reconstructed Y position, uncorrected [cm]',
                  's2_y_3'), np.float32),
                #fifth biggest
                (('fifth S2 time since unix epoch [ns]',
                  's2_time_4'), np.int64),
                (('fifth S2 weighted center time since unix epoch [ns]',
                  's2_center_time_4'), np.int64),
                (('Area of fifth alternate S2 in event [PE]',
                  's2_area_4'), np.float32),]   
        dtype += strax.time_fields
        return dtype
        
    def compute_loop(self, event, peaks):
        result = dict(n_peaks=len(peaks),
                      time=event['time'],
                      endtime=strax.endtime(event))
        if not len(peaks):
            return results
        main_s = dict()
        largest_s = []
        max_number_s1 = 2
        max_number_s2 = 4
        for s_i in [2, 1]:
            if s_i == 2:
                to_store = ['x', 'y']  
            s_mask = peaks['type'] == s_i 
            # For determining the main / alternate S1s, remove all peaks after the main S2 (if there was one). 
            # since these cannot be related to the main S2. This is why S2 finding happened first.
            if s_i == 1 and result[f's2_index'] != -1: 
                s_mask &= peaks['time'] < main_s[2]['time']         
            ss = peaks[s_mask] #peaks of s1
            s_indices = np.arange(len(peaks))[s_mask]

            if not len(ss):
                result[f's{s_i}_index'] = -1
                continue
            largest_i = np.argsort(ss['area'])[::-1]
            result[f's{s_i}_index'] = s_indices[largest_i[0]]
            
            if s_i == 1 :
                if len(ss['area']) < max_number_s1:
                    max_number_s1 = len(ss['area'])
            if s_i == 2 :
                if len(ss['area']) < max_number_s2:
                    max_number_s2 = len(ss['area']) 
            
            irange = range(1,max_number_s1) if s_i == 1 else range(1,max_number_s2)
            for i in irange:
                if len(ss['area']) >= i:
                    index = largest_i[i]
                    x = ss[index]
                    for prop in ['area', 'time', 'center_time'] :
                        result[f's{s_i}_{prop}_{i}'] = x[prop]
                    if s_i == 2:
                        for name in to_store:
                            result[f's2_{name}_{i}'] = x[name]
                else: break
            s = main_s[s_i] = ss[largest_i[0]]
            for prop in ['area', 'time', 'center_time']:
                result[f's{s_i}_{prop}_0'] = s[prop]  
            if s_i == 2:
                for name in to_store:
                    result[f's2_{name}_0'] = s[name]     
                
        return result
