import re
Srbd_amp="ATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGT"
Spbs_amp="TATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTAC"
RdRP_amp="GATGCCACAACTGCTTATGCTAATAGTGTTTTTAACATTTGTCAAGCTGTCACGGCCAATGTTAATGCACTTTTATCTACTGATGGTAACAAAATTGCCGATAAGTATGTCCGCAA"

fileoutSrbd = open("R1R2_Srbd_top_PCR_hits.txt", "w")
fileoutSrbd.write(">Refseq_Srbd_pst_strand" + "\n" + Srbd_amp + "\n")

fileoutSpbs = open("R1R2_Spbs_top_PCR_hits.txt", "w")
fileoutSpbs.write(">Refseq_Spbs_pst_strand" + "\n" + Spbs_amp + "\n")

fileoutRdRP = open("R1R2_RdRP_top_PCR_hits.txt", "w")
fileoutRdRP.write(">Refseq_RdRP_pst_strand" + "\n" + RdRP_amp + "\n")

fileout_combMut = open("important_vars_detailed_V4.txt", "w")
second_info="sample,Srbd_tophit,Spbs_tophit,Srbd_tophit_percentage,Srbd_total_read_count,Spbs_tophit_percentage,Spbs_total_read_count,Srbd_WT_percentage,total_N501Y_percentage,total_E484K_percentage,total_S494P_percentage,Spbs_WT_percentage,total_P681H_percentage,total_P681R_percentage,file"
second_info_l = second_info.split(",")
fileout_combMut.write(second_info + "\n")


with open('SE_fastq_files.txt') as f:
    for line in f:
        filename = line.strip()
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
        rdrp_c_ = 0.001
        Srbd_c_ = 0.001
        Spbs_c_ = 0.001
        N501Yc = 0
        E484Kc =0
        S494Pc =0
        P681Hc =0
        P681Rc =0
        Srbd_WTc =0
        Spbs_WTc =0


        #fileout1 = open((filename[:-4] + '_RdRP_count.txt'), 'w')
        #fileout4 = open((filename[:-4] + '_S_RBD_count.txt'), 'w')
        #fileout8 = open((filename[:-4] + '_S_PBS_count.txt'), 'w')

        long_read_t = long_read_dict.items()
        for key,val in long_read_t :
            if "ACTGCTTATG" in key and "CGATAAGT" in key :
                rdrp_.append((val,key))
                rdrp_c_ = rdrp_c_ + val
            elif "GGTAGCACACCT" in key and "GTTGGTTA" in key:
                Srbd_.append((val, key))
                Srbd_c_ = Srbd_c_+ val
                if "CCCACTTATG" in key: N501Yc=N501Yc+val
                if "TGTTAAAGGT" in key: E484Kc=E484Kc+val
                if "ACAACCATAT" in key: S494Pc=S494Pc+val
                if key in Srbd_amp: Srbd_WTc=Srbd_WTc+val
            elif "TCAGACTCAGAC" in key and "GTGCAGAA" in key:
                Spbs_.append((val, key))
                Spbs_c_ = Spbs_c_+ val
                if "TCTCATCGG" in key: P681Hc=P681Hc+val
                if "TCTCGTCGG" in key: P681Rc=P681Rc+val
                if key in Spbs_amp: Spbs_WTc=Spbs_WTc+val


        comb_mutdic = dict()
        filename_l=filename.split("_")
        comb_mutdic["sample"]=filename_l[1]
        comb_mutdic["file"]=filename
        comb_mutdic["Srbd_WT_percentage"]=("{:.2f}".format(Srbd_WTc/Srbd_c_*100))
        comb_mutdic["Spbs_WT_percentage"]=("{:.2f}".format(Spbs_WTc/Spbs_c_*100))
        comb_mutdic["total_N501Y_percentage"]=("{:.2f}".format(N501Yc/Srbd_c_*100))
        comb_mutdic["total_E484K_percentage"]=("{:.2f}".format(E484Kc/Srbd_c_*100))
        comb_mutdic["total_S494P_percentage"]=("{:.2f}".format(S494Pc/Srbd_c_*100))
        comb_mutdic["total_P681H_percentage"]=("{:.2f}".format(P681Hc/Spbs_c_*100))
        comb_mutdic["total_P681R_percentage"]=("{:.2f}".format(P681Rc/Spbs_c_*100))
        comb_mutdic["Spbs_total_read_count"]=Spbs_c_
        comb_mutdic["Srbd_total_read_count"]=Srbd_c_

        rdrp_.sort(reverse=True)
        rdrp_ltop = rdrp_[:1]
        for val,item in rdrp_ltop :
            if item[:96] not in RdRP_amp :
                if val>2 : fileoutRdRP.write(">" + filename[:18] + "_RdRP_" + str(val) + "of" + str(rdrp_c_) + "\n" + item + "\n")
            #if item in "GATGCCACAACTGCTTATGCTAATAGTGTTTTTAACATTTGTCAAGCTGTCACGGCCAATGTTAATGCACTTTTATCTACTGATGGTAACAAAATTGCCGATAAGTATGTCCGCAA":
                #fileoutMutSum.write(filename + ",WT," + str(val) + ",")
            #else : fileoutMutSum.write(filename + ",new variant," + str(val) + ",")

        #for val,item in rdrp_ :
            #fileout1.write(">" + filename[:18] + "_RdRP_" + str(val) + "of" + str(rdrp_c_) + " " + item + "\n")
        #    fileout1.write(">" + filename[:18] + "_RdRP_" + str("{:.2f}".format(val/rdrp_c_*100)) + "%" + " " + item + "\n")

        Srbd_.sort(reverse=True)
        Srbd_ltop = Srbd_[:1]
        for val,item in Srbd_ltop :
            if item[:87] in Srbd_amp :
                comb_mutdic["Srbd_tophit"]="WT"
                comb_mutdic["Srbd_tophit_percentage"] = ("{:.2f}".format(val/Srbd_c_*100))
            elif item[:87] not in Srbd_amp :
                if val>2 :
                    fileoutSrbd.write(">" + filename[:18] + "_Srbd_" + str(val) + "of" + str(rdrp_c_) + "\n" + item + "\n")
                comb_mutdic["Srbd_tophit"]="-"
                if "CCCACTTATG" in item:
                    comb_mutdic["Srbd_tophit"]=comb_mutdic["Srbd_tophit"] + "N501Y-"
                    comb_mutdic["Srbd_tophit_percentage"] = ("{:.2f}".format(val/Srbd_c_*100))

                if "TGTTAAAGGT" in item:
                    comb_mutdic["Srbd_tophit"]=comb_mutdic["Srbd_tophit"] + "E484K-"
                    comb_mutdic["Srbd_tophit_percentage"] = ("{:.2f}".format(val/Srbd_c_*100))

                if "ACAACCATAT" in item:
                    comb_mutdic["Srbd_tophit"]=comb_mutdic["Srbd_tophit"] + "S494P-"
                    comb_mutdic["Srbd_tophit_percentage"] = ("{:.2f}".format(val/Srbd_c_*100))

                if  "ACAACCATAT" not in item and "TGTTAAAGGT" not in item and "CCCACTTATG" not in item:
                    comb_mutdic["Srbd_tophit"]=comb_mutdic["Srbd_tophit"] + "other variants-"
                    comb_mutdic["Srbd_tophit_percentage"] = ("{:.2f}".format(val/Srbd_c_*100))

        #for val,item in Srbd_ :
            #fileout4.write(">" + filename[:18] + "_Srbd_" + str(val) + "of" + str(rdrp_c_) + " " + item + "\n")

        Spbs_.sort(reverse=True)
        Spbs_ltop = Spbs_[:1]
        for val,item in Spbs_ltop :
            if item not in Spbs_amp :
                if val>2 : fileoutSpbs.write(">" + filename[:18] + "_Spbs_" + str(val) + "of" + str(rdrp_c_) + "\n" + item + "\n")
            if item[:87] in Spbs_amp:
                comb_mutdic["Spbs_tophit"]="WT"
                comb_mutdic["Spbs_tophit_percentage"] = ("{:.2f}".format(val/Spbs_c_*100))
            elif "TCTCATCGG" in item:
                comb_mutdic["Spbs_tophit"]="-P681H-"
                comb_mutdic["Spbs_tophit_percentage"] = ("{:.2f}".format(val/Spbs_c_*100))
            elif "TCTCGTCGG" in item:
                comb_mutdic["Spbs_tophit"]="-P681R-"
                comb_mutdic["Spbs_tophit_percentage"] = ("{:.2f}".format(val/Spbs_c_*100))
            else :
                comb_mutdic["Spbs_tophit"]="-other variants-"
                comb_mutdic["Spbs_tophit_percentage"] = ("{:.2f}".format(val/Spbs_c_*100))
        for items in second_info_l :
            if items in comb_mutdic: fileout_combMut.write(str(comb_mutdic[items]) + ",")
            else: fileout_combMut.write( "NA,")
        fileout_combMut.write( "\n")
        #for val,item in Spbs_ :
            #fileout8.write(">" + filename[:18] + "_Spbs_" + str(val) + "of" + str(rdrp_c_) + " " + item + "\n")
        #fileoutMutSum.write(filename + ",N501Y," + str(N501Yc) +","+ str(Srbd_c_)+ ",Srbd" +"\n")
        #fileoutMutSum.write(filename + ",E484K," + str(E484Kc) +","+ str(Srbd_c_)+ ",Srbd" +"\n")
        #fileoutMutSum.write(filename + ",S494P," + str(S494Pc) +","+ str(Srbd_c_)+ ",Srbd" +"\n")
        #fileoutMutSum.write(filename + ",P681H," + str(P681Hc) +","+ str(Spbs_c_)+ ",Spbs" +"\n")
