fid1 = open("parallel_12/tank_so.log","r")
fid2 = open("serial/tank_so.log","r")

serial = fid2.readlines()
parallel = fid1.readlines()

serialS = 0.
serialE = 0.
id = 0
for line in serial:
    if "Solving over interval" in line and id == 0:
        serialS = float(line[1:9])
        id = 1
    if "Closing Archive" in line:
        serialE = float(line[1:9])
print "serialS,serialE"

parallelS = 0.
parallelE = 0.
id = 0
for line in parallel:
    if "Solving over interval" in line and id == 0:
        parallelS = float(line[1:9])
        id = 1
    if "Closing Archive" in line:
        parallelE = float(line[1:9])

print "Serial start time=" +str(serialS)
print "Serial end time=" +str(serialE)
print "parallel start time=" +str(parallelS)
print "parallel start time=" +str(parallelE)
print "Number of processors=12"
scalEff = (serialE-serialS)/(parallelE-parallelS)/12
print "Scaling efficiency="+str(scalEff)
