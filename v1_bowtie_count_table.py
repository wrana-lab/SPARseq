# 1: open the text file that has the list of all count file names (do not include read2 file names!). Open each count file one by one.
fileout = open("bowtie_count_table_v3.txt", "w")
it=0
fileout.write("samples")
with open('samCountOuts.txt') as f:
    for line in f:
        it=it+1
        countfile = line.strip()
        with open(countfile) as f1:
            for line in f1:
                info = line.strip()
                infol = info.split()
                fileout.write("," + infol[0])
        if it>0 :break
fileout.write("\n")
with open('samCountOuts.txt') as f2:
    for line in f2:
        countfile = line.strip()
        fileout.write(countfile)
        with open(countfile) as f1:
            for line in f1:
                info = line.strip()
                infol = info.split()
                fileout.write("," + infol[1])
        fileout.write("\n")
