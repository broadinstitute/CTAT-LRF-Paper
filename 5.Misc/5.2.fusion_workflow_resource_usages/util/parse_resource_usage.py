#!/usr/bin/env python3

import sys, os, re
import gzip
import glob


def main():


    print("\t".join(["sample", "prog", "exec_time_min", "max_memory_GiB"]))
    
    for monitoring_log_file in glob.glob("*monitoring.log.gz"):

        m = re.match("call-(.*)\\.DepMap_v1v2mrgd_(.*)_monitoring.log.gz", monitoring_log_file)
        if m is not None:
            prog = m.group(1)
            sample_id = m.group(2)

        else:
            raise RuntimeError("no regex match")

        memory = list()

        with gzip.open(monitoring_log_file, 'rt') as fh:
            for line in fh:
                if re.match("\\* Memory usage:", line):
                    pts = line.split()
                    if pts[4] == "GiB":
                        memval = pts[3]
                        memory.append(memval)

        #print(str(memory))
        memory = [float(x) for x in memory]

        max_memory = str(max(memory))

        exec_time = "{:.1f}".format(10 * len(memory) / 60.0)

        print("\t".join([sample_id, prog, exec_time, max_memory]))
    
    
    sys.exit(0)

    



if __name__=='__main__':
    main()
