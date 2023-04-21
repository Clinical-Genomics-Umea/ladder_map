import matplotlib.pyplot as plt

def plot_trace(fsa_obj, trace_name, raw=False):
    if raw:
        plt.plot(fsa_obj.traces_raw[trace_name])
        plt.show()
    else:
        plt.plot(fsa_obj.traces[trace_name])
        plt.show()

