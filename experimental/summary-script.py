"""
This is a completely one-off script.
"""

import os

if __name__ == '__main__':
    outpath = 'out.summary'
    base = '/media/DATAPART2/pileups'
    for filename in os.listdir(base):
        # write the filename
        cmd = 'echo "%s" >> %s' % (filename, outpath)
        print cmd
        os.system(cmd)
        # write the summary
        datapath = os.path.join(base, filename)
        exepath = '/home/argriffi/repos/piledriver/summarize-pileup'
        cmd = 'zcat %s | %s >> %s' % (datapath, exepath, outpath)
        print cmd
        os.system(cmd)

