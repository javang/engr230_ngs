
import csv
import sys
f = open(sys.argv[1])
fo = open(sys.argv[2],"w")

reader =csv.reader(f,delimiter=",")
writer = csv.writer(fo,delimiter=",")
writer.writerow(reader.next())
for row in reader:
    newrow = row
    newrow[4] = float(row[4])**0.5
    writer.writerow(newrow)
f.close()
fo.close()
