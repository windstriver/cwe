# Convert probe post data to csv

probePost = "../postProcessing/probes/0/U"
probeOut = "internal_cell_velocity.csv"

with open(probePost, 'r') as fin:
    with open(probeOut, 'w') as fout:
        for line in fin:
            fout.write(line.replace('(', ' ').replace(')', ' '))
