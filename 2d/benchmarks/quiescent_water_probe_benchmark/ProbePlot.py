import sys
import os.path
from matplotlib import pyplot as plt
import numpy as np

# Proteus native probe
proteus_gauge_file = 'combined_gauge_0_0.5_sample_all.csv'
# Paraview probe
paraview_gauge_file = 'probe_0.5_0.5_0.csv'

global interactive
interactive = True

def parse_proteus_gauge_fields(header):
    def parse_gauge_field(raw_gauge_field):
        gauge_field = raw_gauge_field.strip().rstrip(']')
        field_name, ignore1, x, y, z = gauge_field.split()
        return field_name, (float(x), float(y), float(z))

    assert (header.split(',')[0].strip() == 'time')
    return [parse_gauge_field(gauge_field) for gauge_field in header.split(',')[1:]]


def load_proteus_gauge_data(filename):
    with open(filename) as f:
        proteus_gauge_header = f.readline()
        gauge_fields = parse_proteus_gauge_fields(proteus_gauge_header)
    proteus_timed_gauge_data = np.loadtxt(filename, skiprows=1, delimiter=",")
    return gauge_fields, proteus_timed_gauge_data


def load_paraview_gauge_data(filename):
    return np.loadtxt(filename, skiprows=1, delimiter=",")


def get_theoretical_zero_five(gauge_field, times):
    field, location = gauge_field
    if location != (0.5, 0.5, 0.0):
        raise ValueError("Unexpected location of gauge:", field, "at", location)
    if field == 'p':
        theoretical = np.ones_like(times) * 1.2 * 9.81 * 1.2 + 0.1 * 9.81 * 1000
    elif field in ('u', 'v'):
        theoretical = np.zeros_like(times)
    else:
        raise ValueError("Unrecognize Field", field, "at", location)
    return theoretical


def plot_gauge(times, field, location, measured_data, theoretical_data):
    plt.title(field + " " + str(location))
    plt.xlabel("Time")
    if field == 'p':
        unit = '(Pa)'
        #lim1 = 950
        #lim2 = 1050
    elif field in ('u', 'v'):
        unit = "(m/s)"
        #lim1 = -0.03
        #lim2 = 0.03
    else:
        raise ValueError("Unexpected field, " + field)
    plt.plot(times, measured_data, "b-x", linewidth=2, markersize=3)
    plt.plot(times, theoretical_data, "g-^", linewidth=2, markersize=3)
    plt.grid()
    plt.legend(("Proteus", "Theoretical"))
    plt.ylabel("%s %s" % (field, unit))
    plt.xlim(times[0]-0.1, times[-1]+0.1)
    #plt.ylim(lim1, lim2)
    plt.savefig(str(field) + ".png")
    if interactive:
        plt.show()
    plt.clf()

def main():
    gauge_fields, proteus_timed_gauge_data = load_proteus_gauge_data(proteus_gauge_file)

    # paraview gauge data is currently unplotted
    if os.path.exists(paraview_gauge_file):
        paraview_gauge_data = load_paraview_gauge_data(paraview_gauge_file)
    else:
        paraview_gauge_data = None

    # Split off gauge data
    # First column is time, verified earlier
    times = proteus_timed_gauge_data[:, 0]
    proteus_gauge_data = proteus_timed_gauge_data[:, 1:]
    for gauge_field, measured_gauge_data in zip(gauge_fields, proteus_gauge_data.T):
        field, location = gauge_field
        if location != (0.5, 0.5, 0):
            sys.stderr.write("Ignoring location: " + str(location) + '\n')
            continue
        theoretical_gauge_data = get_theoretical_zero_five(gauge_field, times)
        plot_gauge(times, field, location, measured_gauge_data, theoretical_gauge_data)


if __name__ == '__main__':
    interactive = False
    main()