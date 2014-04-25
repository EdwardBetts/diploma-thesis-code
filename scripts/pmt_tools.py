def pe_pmt(x, y, pedestal, single_pe):
    return (float(x) - pedestal) / (single_pe - pedestal)

def n_mean_pmt(arr, pedestal, single_pe, cutoff=0):
    #cutoff probably not needed
    pe_sum = sum(pe_pmt(i, j, pedestal, single_pe) * j for i, j in enumerate(arr) if i > cutoff)
    events_sum = sum(arr[cutoff:])

    return pe_sum / events_sum
