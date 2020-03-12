import numpy as np

import strax
export, __all__ = strax.exporter()


@export
class DistinctChannels(strax.LoopPlugin):
    __version__ = '0.1.0'
    depends_on = ('event_basics', 'peaks')
    dtype = [
        ('alt_s1_distinct_channels',
         np.int32,
         'Number of PMTs contributing to the secondary S1 '
         'that do not contribute to the main S1'), ] + strax.time_fields

    def compute_loop(self, event, peaks):
        if event['alt_s1_index'] == -1:
            n_distinct = 0
        else:
            s1_a = peaks[event['s1_index']]
            s1_b = peaks[event['alt_s1_index']]
            n_distinct = np.sum(
                (s1_a['area_per_channel'] > 0)
                != (s1_b['area_per_channel'] > 0))

        return dict(alt_s1_distinct_channels=n_distinct,
                    time=event['time'],
                    endtime=event['endtime'])


@export
class EventInfoKr(strax.MergeOnlyPlugin):
    """Alternate version of event_info for Kr and double scatter analysis
      * Uses a different naming convention:
        s1 -> s1_a, alt_s1 -> s1_b, and similarly for s2s;
      * Adds s1_b_distinct_channels, which is relatively expensive to compute
    """

    depends_on = ['event_info', 'distinct_channels']
    provides = 'event_info_kr'

    @staticmethod
    def rename_field(orig_name):
        special_cases = {'alt_cs1': 'cs1_b',
                         'alt_cs2': 'cs2_b',
                         'cs1': 'cs1_a',
                         'cs2': 'cs2_a'}
        if orig_name in special_cases:
            return special_cases[orig_name]

        name = orig_name
        for s_i in [1, 2]:
            if name.startswith(f's{s_i}'):
                name = name.replace(f's{s_i}', f's{s_i}_a')
            if name.startswith(f'alt_s{s_i}'):
                name = name.replace(f'alt_s{s_i}', f's{s_i}_b')
        return name

    def infer_dtype(self):
        self.input_dtype = (
            strax.unpack_dtype(self.deps['event_info'].dtype)
            + [strax.unpack_dtype(self.deps['distinct_channels'].dtype)[0]])
        return [
            ((comment, self.rename_field(name)), dt)
            for (comment, name), dt in self.input_dtype]

    def compute(self, events):
        result = np.zeros(len(events), dtype=self.dtype)
        for (_, name), _ in self.input_dtype:
            result[self.rename_field(name)] = events[name]
        return result
