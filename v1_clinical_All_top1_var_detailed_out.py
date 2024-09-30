fileout_RdRP = open("SparSeq_RdRP_top1_Var_summary_table.txt", "w")
fileout_RdRP.write("sample,well,RdRP_total_read#,RdRP_top1_read#,nt_var,aa_var,present_in_sec_hit,present_in_third_hit,filename" + "\n")
fileout_Srbd = open("SparSeq_Srbd_top1_Var_summary_table.txt", "w")
fileout_Srbd.write("sample,well,Srbd_total_read#,Srbd_top1_read#,nt_var,aa_var,present_in_sec_hit,present_in_third_hit,filename" + "\n")
fileout_Spbs = open("SparSeq_Spbs_top1_Var_summary_table.txt", "w")
fileout_Spbs.write("sample,well,Spbs_total_read#,Spbs_top1_read#,nt_var,aa_var,present_in_sec_hit,present_in_third_hit,filename" + "\n")
#fileout.write("sample info,variation(s) in top-seq, top seq variant(S) present in 2.top-hit?,top seq variant(S) present in 3.top-hit?, top-hit read count, variant(s) in 2.top-hit, 2.top-hit seq variant(S) present in 3.top-hit?, 2.top-hit read count, variants in 3.top-hit, 3.top-hit read count, total amplicon read count" + "\n" )

#fileout1 = open("RdRP_nt_list.txt", "w")
#fileout2 = open("RdRP_aa_list.txt", "w")
#fileout3 = open("Srbd_nt_list.txt", "w")
#fileout4 = open("Srbd_aa_list.txt", "w")
#fileout5 = open("Spbs_nt_list.txt", "w")
#fileout6 = open("Spbs_aa_list.txt", "w")

aa_dic = dict()
with open('codon_aa_table.csv') as f:
    for line in f:
        info = line.strip()
        infol = info.split(",")
        aa_dic[(infol[0])] = infol[1]

RdRP = "GATGCCACAACTGCTTATGCTAATAGTGTTTTTAACATTTGTCAAGCTGTCACGGCCAATGTTAATGCACTTTTATCTACTGATGGTAACAAAATTGCCGATAAGTATGTCCGCAA"
RdRP_nt_dic = dict()
nt_num = 15490
it=0
for RdRP_nt in RdRP :
    RdRP_nt_dic[nt_num]=RdRP[it]
    nt_num=nt_num+1
    it=it+1
RdRP_aa_dic = dict()
pos=0
aa_num=684 #adjusted from 575 to 684 to match ref
for RdRP_aa in RdRP :
    if pos+3<len(RdRP):
        codon = RdRP[pos:(pos+3)]
        RdRP_aa_dic[aa_num] = aa_dic[codon]
        pos = pos+3
        aa_num = aa_num+1
    else : break
#for val in RdRP_nt_dic : fileout1.write(str(val) + "," + RdRP_nt_dic[val] + "\n")
#for val in RdRP_aa_dic : fileout2.write(str(val) + "," + RdRP_aa_dic[val]+ "\n")

Srbd = "ATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGT"
Srbd_nt_dic = dict()
nt_num = 22980
it=0
for Srbd_nt in Srbd :
    Srbd_nt_dic[nt_num]=Srbd[it]
    nt_num=nt_num+1
    it=it+1
Srbd_aa_dic = dict()
pos=2
aa_num=474
for Srbd_aa in Srbd :
    if pos+3<len(Srbd):
        codon = Srbd[pos:(pos+3)]
        Srbd_aa_dic[aa_num] = aa_dic[codon]
        pos = pos+3
        aa_num = aa_num+1
    else : break
#for val in Srbd_nt_dic : fileout3.write(str(val) + "," + Srbd_nt_dic[val]+ "\n")
#for val in Srbd_aa_dic : fileout4.write(str(val) + "," + Srbd_aa_dic[val]+ "\n")

Spbs = "TATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTAC"
Spbs_nt_dic = dict()
nt_num = 23571
it=0
for Spbs_nt in Spbs :
    Spbs_nt_dic[nt_num]=Spbs[it]
    nt_num=nt_num+1
    it=it+1
Spbs_aa_dic = dict()
pos=2
aa_num=671
for Spbs_aa in Spbs :
    if pos+3<len(Spbs):
        codon = Spbs[pos:(pos+3)]
        Spbs_aa_dic[aa_num] = aa_dic[codon]
        pos = pos+3
        aa_num = aa_num+1
    else : break
#for val in Spbs_nt_dic : fileout5.write(str(val) + "," + Spbs_nt_dic[val]+ "\n")
#for val in Spbs_aa_dic : fileout6.write(str(val) + "," + Spbs_aa_dic[val]+ "\n")

RdRP_var_dic = dict()
Srbd_var_dic = dict()
Spbs_var_dic = dict()

with open('SE_fastq_files.txt') as f:
    for line in f:
        filename = line.strip()
        filename_l = filename.split("_")
        #Part-1 binning sequences and creating "sequence - readcount" file
        #short_read_dict = dict()
        long_read_dict = dict()
        it=2
        with open(filename) as f:
            for line in f:
                it=it+1
                if it%4==0 :
                    rseq = line.strip()
                    adap_pos = rseq.find("AGATCGGAAG")
                    seq = rseq[:adap_pos]
                    if len(seq)>90 :
                        if seq not in long_read_dict : long_read_dict[seq] = 1
                        else : long_read_dict[seq] = 1+long_read_dict[seq]
                    #else:
                        #if seq not in short_read_dict : short_read_dict[seq] = 1
                        #else : short_read_dict[seq] = 1+short_read_dict[seq]
        #Part-2 binning sequences with read counts into PCR pools and sorting writing top
        rdrp_ = list()
        Srbd_ = list()
        Spbs_ = list()
        rdrp_c_ = 0
        Srbd_c_ = 0
        Spbs_c_ = 0
        #fileout1 = open((filename[:-4] + '_RdRP_count.txt'), 'w')
        #fileout4 = open((filename[:-4] + '_S_RBD_count.txt'), 'w')
        #fileout8 = open((filename[:-4] + '_S_PBS_count.txt'), 'w')

        long_read_t = long_read_dict.items()
        for key,val in long_read_t :
            if "ACTGCTTATG" in key and "CGATAAGT" in key:
                rdrp_.append((val,key))
                rdrp_c_ = rdrp_c_ + val
            elif "GGTAGCACACCT" in key and "GTTGGTTA" in key:
                Srbd_.append((val, key))
                Srbd_c_ = Srbd_c_+ val
            elif "TCAGACTCAGAC" in key and "GTGCAGAA" in key:
                Spbs_.append((val, key))
                Spbs_c_ = Spbs_c_+ val

        rdrp_.sort(reverse=True)
        var_cond=0
        for val,item in rdrp_[:1] :
            if val>3 :
                seq=item[2:96]
                nt_num = 15492
                for nt in seq :
                    if nt_num<15587 and RdRP_nt_dic[nt_num]!= nt :
                        fileout_RdRP.write(filename_l[1] + "," + filename_l[4] +","+ str(rdrp_c_) + "," + str(val) + "," + RdRP_nt_dic[nt_num]+str(nt_num)+nt+",")
                        var_aa_num = 684+int((nt_num-15490)/3) #adjust 575 to 684 to match ref
                        var_aa_start = ((var_aa_num-684)*3)
                        var_aa_end = var_aa_start+3
                        ncodon = seq[var_aa_start:var_aa_end]
                        if "N" not in ncodon:fileout_RdRP.write(RdRP_aa_dic[var_aa_num]+str(var_aa_num)+aa_dic[ncodon]+",")
                        else: fileout_RdRP.write(RdRP_aa_dic[var_aa_num]+str(var_aa_num)+"NA"+",")
                        for val2,item2 in rdrp_[1:2] :
                            seq2=item2[2:96]
                            if ncodon == seq2[var_aa_start:var_aa_end]:
                                fileout_RdRP.write("yes#" + str(val2) + ",")
                            else: fileout_RdRP.write("no#" + str(val2) + ",")
                        for val3,item3 in rdrp_[2:3] :
                            seq3=item3[2:96]
                            if ncodon == seq3[var_aa_start:var_aa_end]:
                                fileout_RdRP.write("yes#" + str(val3) + "," + filename + "\n")
                            else: fileout_RdRP.write("no#" + str(val3) + "," + filename + "\n")
                        var_cond =var_cond+1
                    nt_num=nt_num+1
                if var_cond==0 : fileout_RdRP.write(filename_l[1] + "," + filename_l[4] +","+ str(rdrp_c_) + "," + str(val) + ",WT,WT,NA,NA," + filename + "\n")
            else: fileout_RdRP.write(filename_l[1] + "," + filename_l[4] +","+ str(rdrp_c_) + "," + str(val) + ",low,low,low,low," + filename + "\n")
        Srbd_.sort(reverse=True)
        for val,item in Srbd_[:1] :
            if val>3 :
                seq=item[2:87]
                nt_num = 22982
                var_cond=0
                for nt in seq :
                    if nt_num<23068 and Srbd_nt_dic[nt_num]!= nt :
                        fileout_Srbd.write(filename_l[1] + "," + filename_l[4] +","+ str(Srbd_c_) + "," + str(val) + "," + Srbd_nt_dic[nt_num]+str(nt_num)+nt+",")
                        var_aa_num = 474+int((nt_num-22982)/3)
                        var_aa_start = ((var_aa_num-474)*3)
                        var_aa_end = var_aa_start+3
                        ncodon = seq[var_aa_start:var_aa_end]
                        if "N" not in ncodon:fileout_Srbd.write(Srbd_aa_dic[var_aa_num]+str(var_aa_num)+aa_dic[ncodon]+",")
                        else:fileout_Srbd.write(Srbd_aa_dic[var_aa_num]+str(var_aa_num)+"NA"+",")
                        for val2,item2 in Srbd_[1:2] :
                            seq2=item2[2:87]
                            if ncodon == seq2[var_aa_start:var_aa_end]:
                                fileout_Srbd.write("yes#" + str(val2) + ",")
                            else: fileout_Srbd.write("no#" + str(val2) + ",")
                        for val3,item3 in Srbd_[2:3] :
                            seq3=item3[2:87]
                            if ncodon == seq3[var_aa_start:var_aa_end]:
                                fileout_Srbd.write("yes#" + str(val3) + "," + filename + "\n")
                            else: fileout_Srbd.write("no#" + str(val3) + "," + filename + "\n")
                        var_cond =var_cond+1
                    nt_num=nt_num+1
                if var_cond==0 : fileout_Srbd.write(filename_l[1] + "," + filename_l[4] +","+ str(Srbd_c_) + "," + str(val) + ",WT,WT,NA,NA," + filename + "\n")
            else: fileout_Srbd.write(filename_l[1] + "," + filename_l[4] +","+ str(Srbd_c_) + "," + str(val) + ",low,low,low,low," + filename + "\n")

        Spbs_.sort(reverse=True)
        for val,item in Spbs_[:1] :
            if val>3 :
                seq=item[2:87]
                nt_num = 23573
                var_cond=0
                for nt in seq :
                    if nt_num<23662 and Spbs_nt_dic[nt_num]!= nt :
                        fileout_Spbs.write(filename_l[1] + "," + filename_l[4] + ","+ str(Spbs_c_) + "," + str(val) + ","+Spbs_nt_dic[nt_num]+str(nt_num)+nt+",")
                        var_aa_num = 671+int((nt_num-23573)/3)
                        var_aa_start = ((var_aa_num-671)*3)
                        var_aa_end = var_aa_start+3
                        ncodon = seq[var_aa_start:var_aa_end]
                        if "N" not in ncodon:fileout_Spbs.write(Spbs_aa_dic[var_aa_num]+str(var_aa_num)+aa_dic[ncodon]+",")
                        else:fileout_Spbs.write(Spbs_aa_dic[var_aa_num]+str(var_aa_num)+"NA"+",")
                        for val2,item2 in Spbs_[1:2] :
                            seq2=item2[2:87]
                            if ncodon == seq2[var_aa_start:var_aa_end]:
                                fileout_Spbs.write("yes#" + str(val2) + ",")
                            else: fileout_Spbs.write("no#" + str(val2) + ",")
                        for val3,item3 in Spbs_[2:3] :
                            seq3=item3[2:87]
                            if ncodon == seq3[var_aa_start:var_aa_end]:
                                fileout_Spbs.write("yes#" + str(val3) + "," + filename + "\n")
                            else: fileout_Spbs.write("no#" + str(val3) + "," + filename + "\n")
                        var_cond =var_cond+1
                    nt_num=nt_num+1
                if var_cond==0 : fileout_Spbs.write(filename_l[1] + "," + filename_l[4] +","+ str(Spbs_c_) + "," + str(val) + ",WT,WT,NA,NA," + filename + "\n")
            else: fileout_Spbs.write(filename_l[1] + "," + filename_l[4] +","+ str(Spbs_c_) + "," + str(val) + ",low,low,low,low," + filename + "\n")
