import pysam
import argparse
from pyfaidx import Fasta
import time
import copy
from PIL import Image
import numpy as np
import os

#the function to expand cigar and record mismatch
def modify_cigar(cigar,reverse,pos,read_str):
    cigar_str=''
    i=0
    S_clipping=[0,0]
    sign=0
    #this step to expand cigar and record S length
    while i<len(cigar):
        num=''
        while cigar[i].isdigit():
            num=num+cigar[i]
            i+=1
        if sign==0 and cigar[i]=='S':
            S_clipping[0]=int(num)
        if sign==1 and cigar[i]=='S':
            S_clipping[1] = int(num)
        sign=1
        if (cigar[i]!='S' and cigar[i]!='H'):
            cigar_str=cigar_str+cigar[i]*int(num)
        i+=1
    #this step is to get sequence of reference and read
    ref=genes[chr_name][pos:pos+len(cigar_str)]
    read=read_str[S_clipping[0]:len(read_str)-S_clipping[1]]
    cigar_str2=''
    pos=0
    pos1=0
    pos2=0
    while pos<len(cigar_str) and pos1<len(read) and pos2<len(ref):
        if cigar_str[pos]=='M':
            if str(ref[pos2]).upper()==read[pos1].upper():
                cigar_str2+='M'
            else:
                cigar_str2 += 'X'
            pos1 += 1
            pos2 += 1
            pos += 1
        elif cigar_str[pos]=='D':
            cigar_str2 += 'D'
            pos2+=1
            pos += 1
        elif cigar_str[pos]=='I':
            num=0
            while pos1<len(read) and pos<len(cigar_str) and cigar_str[pos]=='I':
                pos1 += 1
                pos += 1
                num+=1
            #INS don't take over a pos in cigar, I delete some m and insert some I. 50I will replace 1 m
            if num>=10:
                I_NUM = int(num / num_i_to_m) + 1
                cigar_str2 = cigar_str2[:-I_NUM]
                cigar_str2 = cigar_str2 + 'I' * I_NUM
        else:
            pos1+=1
            pos += 1
    #this step is to record whether the alignment is reverse
    if reverse==False:
        return cigar_str2.upper()
    else:
        return cigar_str2.lower()

#the function to merge cigar of same read for the read have over 2 alignments
def merge_cigar(cigar1, cigar2, pos1, pos2):
    if pos1<pos2:
        n=pos2-pos1-len(cigar1)
        cigar_str=cigar1+'0'*n+cigar2
        res=[pos1,cigar_str]
    else:
        n = pos1 - pos2 - len(cigar2)
        cigar_str = cigar2 + '0' * n + cigar1
        res = [pos2, cigar_str]
    return res

#the function to read bam file and record cigar
def read_cigarfile(samfile,genes,chr_name,pos1,pos2):
    read_cigar = {}
    read_num=0
    for r in samfile.fetch(chr_name,pos1-500,pos2+500):
        if r.mapq > 20 :
            read_num+=1
            read_name = r.query_name
            if read_name in read_cigar:
                pos1 = read_cigar[read_name][0]
                cigar1 = read_cigar[read_name][1]
                pos2=r.reference_start
                cigar2 = modify_cigar(r.cigarstring, r.is_reverse, r.reference_start, r.query_sequence)
                read_cigar[read_name] = merge_cigar(cigar1, cigar2, pos1, pos2)
            else:
                read_cigar[read_name] = [r.reference_start,
                                         modify_cigar(r.cigarstring, r.is_reverse, r.reference_start, r.query_sequence)]
            if read_num>=200:
                break
    return read_cigar

#the function to transfer cigar into matrix
def create_matrix(read_cigar,total_read):
    matrix = np.zeros((len(read_cigar), matrix_len))
    i=0
    for key in sorted(total_read):
        read_s = int(read_cigar[key][0])
        read_e = int(read_cigar[key][0]) + len(read_cigar[key][1])
        # to make sure whether matrix should start with pos 0
        matrix_s = max(0, read_s - WIN_S)
        matrix_e = min(matrix_len, read_e - WIN_S)
        for pos in range(0, matrix_e - matrix_s):
            if WIN_S - int(read_cigar[key][0]) > 0:
                clipping = WIN_S - int(read_cigar[key][0])
            else:
                clipping = 0
            #D is 2, X is 4, 0 is null
            if read_cigar[key][1][clipping + pos] == 'D' or read_cigar[key][1][clipping + pos] == 'd':
                matrix[i][matrix_s + pos] = 2
            elif read_cigar[key][1][clipping + pos] == 'X' or read_cigar[key][1][clipping + pos] == 'x':
                matrix[i][matrix_s + pos] = 4
            elif read_cigar[key][1][clipping + pos] == 'I' or read_cigar[key][1][clipping + pos] == 'i':
                matrix[i][matrix_s + pos] = 5
            elif read_cigar[key][1][clipping + pos].isupper():
                matrix[i][matrix_s + pos] = 1
            elif read_cigar[key][1][clipping + pos] == '0':
                matrix[i][matrix_s + pos] = 0
            else:
                matrix[i][matrix_s + pos] = 3
        i += 1
    return matrix

#the function to create a image
def create_image(matrix,Image_name,n):
    c = Image.new("RGB", (n*matrix.shape[0], Image_len))
    image_s = int(matrix_len / 2) - int(Image_len / 2)
    image_e = int(matrix_len / 2) + int(Image_len / 2)
    for i in range(0, matrix.shape[0]):
        for j in range(image_s, image_e):
            for k in range(0,n):
                # D is 2, X is 4, 0 is null
                if matrix[i, j] == 0:
                    c.putpixel([i*n+k, j - image_s], (0, 0, 0))
                elif matrix[i, j] == 2:
                    c.putpixel([i*n+k, j - image_s], (0, 0, 0))
                elif matrix[i, j] == 3:
                    c.putpixel([i*n+k, j - image_s], (0,200,255))
                elif matrix[i, j] == 4:
                    c.putpixel([i*n+k, j - image_s], (0,255,0))
                elif matrix[i, j] == 5:
                    c.putpixel([i*n+k, j - image_s], (255,0,0))
                else:
                    c.putpixel([i*n+k, j - image_s], (255, 255, 0))
    c.save('./SV_image/'+Image_name)

def check_ref(ref):
    if len(ref)<=60 and len(ref)>=30:
        A_num = 0
        G_num = 0
        C_num = 0
        T_num = 0
        for base in ref:
            if base == 'a' or base == 'A':
                A_num += 1
            elif base == 'g' or base == 'G':
                G_num += 1
            elif base == 'c' or base == 'C':
                C_num += 1
            elif base == 't' or base == 'T':
                T_num += 1
            A_rate = float(A_num) / float(len(ref))
            G_rate = float(G_num) / float(len(ref))
            C_rate = float(C_num) / float(len(ref))
            T_rate = float(T_num) / float(len(ref))
        if A_rate > 0.9 or G_rate > 0.9 or C_rate > 0.9 or T_rate > 0.9:
            return False
        else:
            return True
    else:
        return True


def calculate_SV_num(read):
    SV_num=0
    for SV in read:
        if SV[4]=='DEL' or SV[4]=='INS':
            SV_num+=1
        elif SV[4]=='INV':
            SV_num-=1
    return SV_num

def calculate_read_len(cigar):
    read_len=0
    i = 0
    while i < len(cigar):
        num = ''
        while cigar[i].isdigit():
            num = num + cigar[i]
            i += 1
        num = int(num)
        if cigar[i] == 'M' or cigar[i] == 'I' or cigar[i] == 'S':
            read_len+=num
        i += 1
    return read_len

def calculate_ref_len(cigar):
    ref_len=0
    i = 0
    while i < len(cigar):
        num = ''
        while cigar[i].isdigit():
            num = num + cigar[i]
            i += 1
        num = int(num)
        if cigar[i] == 'M' or cigar[i] == 'D':
            ref_len+=num
        i += 1
    return ref_len

parser = argparse.ArgumentParser()
parser.add_argument("bam", help='bam file')
parser.add_argument("fasta", help='fasta file')
args = parser.parse_args()
localtime = time.asctime(time.localtime(time.time()))
print ('Start',localtime)
merge_distance=1000
low_rate=0.8
high_rate=1.2
mini_SV_len=10
coverage=[10]
genes = Fasta(args.fasta)

main_chr=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
          'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrY','chrX',
          '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12','13', '14', '15', '16', '17',
          '18', '19', '20', '21', '22', 'Y', 'X']

FindCSV=[]



samfile=pysam.AlignmentFile(args.bam, 'rb')
#Detect the potential SV
record=[]
for r in samfile.fetch():
  if r.mapq>20 and r.is_supplementary==False:
      chr = r.reference_name
      start_pos=r.reference_start
      cigar=r.cigarstring
      read_str=r.query_sequence
      i=0
      pos=start_pos
      pos1=0
      pos2=0
      while i<len(cigar):
          num=''
          while cigar[i].isdigit():
              num=num+cigar[i]
              i+=1
          num=int(num)
          if cigar[i]=='M':
              pos=pos+num
          elif cigar[i]=='D':
              pos1=pos
              pos=pos+num
              pos2=pos
              if num>=30:
                 record.append([chr, pos1, pos2])
          elif cigar[i] == 'I':
              if num >= 30:
                 record.append([chr, pos, pos])
          i+=1

      sign = ""
      for tag in r.tags:
          if str(tag[0]) == "SA":
              sign = "SA"
              supple_align = tag
      # judge the two alignment
      if sign == "SA":
          # check primary alignment
          chr1 = r.reference_name
          cigar1 = r.cigarstring
          read = r.query_sequence
          other_align = supple_align[1].split(";")
          supple = other_align[0].split(",")
          i = 1
          while i < (len(other_align) - 1):
              m = other_align[i].split(",")
              if int(m[4]) > int(supple[4]):
                  supple = m
              i += 1
          chr2 = supple[0]
          # judge the primary alignment and supplementary alignment pos
          if r.reference_start <= int(supple[1]) and chr1 == chr2:
              pos1 = r.reference_end
              pos2 = int(supple[1])
              direct2 = supple[2]
              cigar2 = supple[3]
              s2 = cigar2.split("S")[0]
              if 'M' in s2:
                  s2 = 0
              else:
                  s2 = int(s2)
              if r.is_reverse == True:
                  direct1 = '-'
              else:
                  direct1 = '+'
              read_pos1 = r.query_alignment_end
              read_pos2 = s2
              DEL_len = (read_pos1 - read_pos2) - (pos1 - pos2)
              if (DEL_len >= 200 or abs(pos1-pos2)>=200) and direct2 == direct1:
                  record.append([chr1, min(pos1, pos2), max(pos1, pos2)])

          elif chr1 == chr2:
              pos1 = r.reference_start
              pos2 = int(supple[1])
              direct2 = supple[2]
              cigar2 = supple[3]
              r2_len = 0
              i = 0
              while i < len(cigar2):
                  num = ''
                  while cigar2[i].isdigit():
                      num = num + cigar2[i]
                      i += 1
                  num = int(num)
                  if cigar2[i] == 'M' or cigar2[i] == 'D':
                      r2_len = r2_len + num
                  i += 1
              s2 = 0
              if cigar2[-1] == 'S' or cigar2[-1] == 'H':
                  s2 = num
              pos2 = pos2 + r2_len
              if r.is_reverse == True:
                  direct1 = '-'
              else:
                  direct1 = '+'
              read_pos2 = len(read) - s2
              read_pos1 = r.query_alignment_start
              DEL_len = (read_pos2 - read_pos1) - (pos2 - pos1)
              if (DEL_len >= 200 or abs(pos1-pos2)>=200) and direct2 == direct1:
                  record.append([chr1, min(pos1, pos2), max(pos1, pos2)])
samfile.close()
#merge the candidate SV and identify the candidate region
record2=[]
for site in record:
    if site[2]-site[1]<=5000:
        record2.append(site)
        #w3.write(str(site[0])+' '+str(site[1])+' '+str(site[2])+'\n')

candidate_region=[]
intervals =list(sorted(record2))
chr=intervals[0][0]
low = intervals[0][1]
high = intervals[0][2]
read_num=1
for i in range(1, len(intervals)):
    new_region_len=intervals[i][2]-intervals[i][1]
    if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr:
        low = min(intervals[i][1],low)
        high = max(intervals[i][2],high)
        read_num += 1
    else:
        candidate_region.append([chr,low, high,read_num])
        chr=intervals[i][0]
        low = intervals[i][1]
        high = intervals[i][2]
        read_num=1
candidate_region.append([chr,low, high,read_num])
newres=[]
for site in candidate_region:
    if int(site[3])>=3:
        newres.append([site[0],site[1],site[2]])
        #w4.write(str(site[0]) + ' ' + str(site[1]) + ' ' + str(site[2]) + '\n')
SV_region=newres
#detect SV in candidate SV region
newSV=[]
samfile=pysam.AlignmentFile(args.bam, 'rb')
for region in SV_region:
    chr = region[0]
    sv_pos1 = region[1]
    sv_pos2 = region[2]
    range1 = sv_pos1 - 500
    range2 = sv_pos2 + 500
    read_num = 0
    read_SV_len = []
    for r in samfile.fetch(chr, sv_pos1, sv_pos2):
        read_name = r.query_name
        if read_num>=150:
            break
        if r.mapq > 20 and r.is_supplementary == False and r.reference_start < range1 and r.reference_end > range2:
            read_num += 1
            read_SV = []
            DEL_NUM = 0
            INS_NUM = 0
            DEL_len = 0
            INS_len = 0
            chr = r.reference_name
            cigar = r.cigarstring
            # read_str = r.query_sequence
            read_name = r.query_name
            i = 0
            pos = r.reference_start
            pos1 = 0
            pos2 = 0
            while i < len(cigar):
                num = ''
                while cigar[i].isdigit():
                    num = num + cigar[i]
                    i += 1
                num = int(num)
                if cigar[i] == 'M':
                    pos = pos + num
                elif cigar[i] == 'D':
                    pos1 = pos
                    pos = pos + num
                    pos2 = pos
                    if num >= mini_SV_len and max(pos1, range1) < min(pos2, range2):
                        read_SV.append([chr, pos1, pos2, num, 'DEL'])
                        DEL_len += num
                        DEL_NUM += 1
                elif cigar[i] == 'I':
                    if num >= mini_SV_len and pos > range1 and pos < range2:
                        INS_len += num
                        read_SV.append([chr, pos, pos, num, 'INS'])
                        INS_NUM += 1
                i += 1
            read_SV_len.append([DEL_len - INS_len, DEL_len, DEL_NUM, INS_len, INS_NUM, read_SV, read_name])

        elif r.mapq > 20 and r.is_supplementary == False:
            read_num+=1
            read_SV = []
            DEL_NUM = 0
            INS_NUM = 0
            DEL_len = 0
            INS_len = 0
            read_name = r.query_name
            sign = ""
            for tag in r.tags:
                if str(tag[0]) == "SA":
                    sign = "SA"
                    supple_align = tag
            # judge the two alignment
            if sign == "SA":
                # check primary alignment
                chr1 = r.reference_name
                cigar1 = r.cigarstring
                read = r.query_sequence
                other_align = supple_align[1].split(";")
                supple = other_align[0].split(",")
                #supple is another alignment, supple2 is middle part of INV
                if r.is_reverse == True:
                    direct1 = '-'
                else:
                    direct1 = '+'
                supple2 = None
                MAPQ_value = 0
                i = 1
                while i < (len(other_align) - 1):
                    m = other_align[i].split(",")
                    if calculate_read_len(m[3]) > calculate_read_len(supple[3]) and chr1==m[0]:
                        supple = m
                    if m[2] != direct1 and int(m[4]) > MAPQ_value and chr1==m[0]:
                        supple2 = m
                        MAPQ_value = int(m[4])
                    i += 1

                chr2 = supple[0]
                if chr1!=chr2:
                    continue
                supple_end=int(supple[1])+calculate_ref_len(supple[3])

                if supple2:
                    INV_len=calculate_ref_len(supple2[3])
                    supple2_end = int(supple2[1]) + INV_len
                # judge the primary alignment and supplementary alignment pos
                if r.reference_start <= int(supple[1]) and chr1 == chr2 and\
                        r.reference_start<=range1 and supple_end>=range2:
                    if supple2 and int(supple2[1])>=range1 and supple2_end<=range2:
                        read_SV.append([chr, int(supple2[1]), supple2_end, INV_len, 'INV'])
                    cigar = r.cigarstring
                    i = 0
                    pos = r.reference_start
                    pos1 = 0
                    pos2 = 0
                    while i < len(cigar):
                        num = ''
                        while cigar[i].isdigit():
                            num = num + cigar[i]
                            i += 1
                        num = int(num)
                        if cigar[i] == 'M':
                            pos = pos + num
                        elif cigar[i] == 'D':
                            pos1 = pos
                            pos = pos + num
                            pos2 = pos
                            if num >= mini_SV_len and max(pos1, range1) < min(pos2, range2):
                                read_SV.append([chr, pos1, pos2, num, 'DEL'])
                                DEL_len += num
                                DEL_NUM += 1
                        elif cigar[i] == 'I':
                            if num >= mini_SV_len and pos > range1 and pos < range2:
                                INS_len += num
                                read_SV.append([chr, pos, pos, num, 'INS'])
                                INS_NUM += 1
                        i += 1

                    pos1 = r.reference_end
                    pos2 = int(supple[1])
                    direct2 = supple[2]
                    cigar2 = supple[3]
                    s2 = cigar2.split("S")[0]
                    if 'M' in s2:
                        s2 = 0
                    else:
                        s2 = int(s2)
                    if r.is_reverse == True:
                        direct1 = '-'
                    else:
                        direct1 = '+'
                    read_pos1 = r.query_alignment_end
                    read_pos2 = s2
                    SV_len = (read_pos1 - read_pos2) - (pos1 - pos2)
                    ref_dis= pos2 - pos1
                    read_dis= read_pos2 - read_pos1
                    if ref_dis>=50:
                        read_SV.append([chr, pos1, pos2, ref_dis, 'DEL'])
                        DEL_len = DEL_len+ref_dis
                        DEL_NUM += 1
                        if read_dis>=50:
                            read_SV.append([chr, pos1, pos2, read_dis, 'INS'])
                            INS_len = INS_len+read_dis
                            INS_NUM += 1
                    else:
                        if SV_len > 200 and min(pos1, pos2) > range1 and max(pos1, pos2) < range2:
                            read_SV.append([chr, pos1, pos2, SV_len, 'DEL'])
                            DEL_len = SV_len
                            DEL_NUM += 1
                        elif SV_len < -200 and min(pos1, pos2) > range1 and max(pos1, pos2) < range2:
                            read_SV.append([chr, pos1, pos2, abs(SV_len), 'INS'])
                            INS_len = abs(SV_len)
                            INS_NUM += 1

                    cigar = supple[3]
                    i = 0
                    pos = int(supple[1])
                    pos1 = 0
                    pos2 = 0
                    while i < len(cigar):
                        num = ''
                        while cigar[i].isdigit():
                            num = num + cigar[i]
                            i += 1
                        num = int(num)
                        if cigar[i] == 'M':
                            pos = pos + num
                        elif cigar[i] == 'D':
                            pos1 = pos
                            pos = pos + num
                            pos2 = pos
                            if num >= mini_SV_len and max(pos1, range1) < min(pos2, range2):
                                read_SV.append([chr, pos1, pos2, num, 'DEL'])
                                DEL_len += num
                                DEL_NUM += 1
                        elif cigar[i] == 'I':
                            if num >= mini_SV_len and pos > range1 and pos < range2:
                                INS_len += num
                                read_SV.append([chr, pos, pos, num, 'INS'])
                                INS_NUM += 1
                        i += 1
                    read_SV_len.append([DEL_len - INS_len, DEL_len, DEL_NUM, INS_len, INS_NUM, read_SV, read_name])

                elif chr1 == chr2 and int(supple[1])<=range1 and r.reference_end>=range2:
                    cigar = supple[3]
                    if supple2 and int(supple2[1])>=range1 and supple2_end<=range2:
                        read_SV.append([chr, int(supple2[1]), supple2_end, INV_len, 'INV'])
                    i = 0
                    pos = int(supple[1])
                    pos1 = 0
                    pos2 = 0
                    while i < len(cigar):
                        num = ''
                        while cigar[i].isdigit():
                            num = num + cigar[i]
                            i += 1
                        num = int(num)
                        if cigar[i] == 'M':
                            pos = pos + num
                        elif cigar[i] == 'D':
                            pos1 = pos
                            pos = pos + num
                            pos2 = pos
                            if num >= mini_SV_len and max(pos1, range1) < min(pos2, range2):
                                read_SV.append([chr, pos1, pos2, num, 'DEL'])
                                DEL_len += num
                                DEL_NUM += 1
                        elif cigar[i] == 'I':
                            if num >= mini_SV_len and pos > range1 and pos < range2:
                                INS_len += num
                                read_SV.append([chr, pos, pos, num, 'INS'])
                                INS_NUM += 1
                        i += 1

                    pos1 = r.reference_start
                    pos2 = int(supple[1])
                    direct2 = supple[2]
                    cigar2 = supple[3]
                    r2_len = 0
                    i = 0
                    while i < len(cigar2):
                        num = ''
                        while cigar2[i].isdigit():
                            num = num + cigar2[i]
                            i += 1
                        num = int(num)
                        if cigar2[i] == 'M' or cigar2[i] == 'D':
                            r2_len = r2_len + num
                        i += 1
                    s2 = 0
                    if cigar2[-1] == 'S' or cigar2[-1] == 'H':
                        s2 = num
                    pos2 = pos2 + r2_len
                    if r.is_reverse == True:
                        direct1 = '-'
                    else:
                        direct1 = '+'
                    read_pos2 = len(read) - s2
                    read_pos1 = r.query_alignment_start
                    SV_len = (read_pos2 - read_pos1) - (pos2 - pos1)
                    ref_dis = pos1 - pos2
                    read_dis = read_pos1 - read_pos2
                    if ref_dis>=50:
                        read_SV.append([chr, pos2, pos1, ref_dis, 'DEL'])
                        DEL_len = DEL_len+ref_dis
                        DEL_NUM += 1
                        if read_dis>=50:
                            read_SV.append([chr, pos2, pos1, read_dis, 'INS'])
                            INS_len = INS_len+read_dis
                            INS_NUM += 1
                    else:
                        if SV_len > 200 and min(pos1, pos2) > range1 and max(pos1, pos2) < range2:
                            read_SV.append([chr, pos2, pos1, SV_len, 'DEL'])
                            DEL_len = SV_len
                            DEL_NUM += 1
                        elif SV_len < -200 and min(pos1, pos2) > range1 and max(pos1, pos2) < range2:
                            read_SV.append([chr, pos2, pos1, abs(SV_len), 'INS'])
                            INS_len = abs(SV_len)
                            INS_NUM += 1

                    cigar = r.cigarstring
                    i = 0
                    pos = r.reference_start
                    pos1 = 0
                    pos2 = 0
                    while i < len(cigar):
                        num = ''
                        while cigar[i].isdigit():
                            num = num + cigar[i]
                            i += 1
                        num = int(num)
                        if cigar[i] == 'M':
                            pos = pos + num
                        elif cigar[i] == 'D':
                            pos1 = pos
                            pos = pos + num
                            pos2 = pos
                            if num >= mini_SV_len and max(pos1, range1) < min(pos2, range2):
                                read_SV.append([chr, pos1, pos2, num, 'DEL'])
                                DEL_len += num
                                DEL_NUM += 1
                        elif cigar[i] == 'I':
                            if num >= mini_SV_len and pos > range1 and pos < range2:
                                INS_len += num
                                read_SV.append([chr, pos, pos, num, 'INS'])
                                INS_NUM += 1
                        i += 1
                    read_SV_len.append([DEL_len - INS_len, DEL_len, DEL_NUM, INS_len, INS_NUM, read_SV, read_name])

    res_len = []
    if len(read_SV_len) <= 2:
        continue
    intervals = list(sorted(read_SV_len))
    jizhi=int(max(len(intervals)/20,1))
    intervals=intervals[jizhi:-jizhi]
    read_cluster=[]
    for every_read in intervals:
        read_cluster.append([every_read[0],1,[[every_read[0],every_read[5]]]])
    #chm13 is 1
    while len(read_cluster)>2:
        read_dis=[]
        i=1
        while i<len(read_cluster):
            read_dis.append(read_cluster[i][0] - read_cluster[i - 1][0])
            i += 1
        mini_sign=0
        mini_dis=read_dis[0]
        j=1
        while j < len(read_dis):
            if read_dis[j] < mini_dis:
                mini_dis = read_dis[j]
                mini_sign = j
            j += 1

        cluster_read_num = read_cluster[mini_sign][1] + read_cluster[mini_sign + 1][1]
        cluster_len = (read_cluster[mini_sign][0] * read_cluster[mini_sign][1] + \
                       read_cluster[mini_sign + 1][0] * read_cluster[mini_sign + 1][1]) / cluster_read_num
        cluster_SV = read_cluster[mini_sign][2] + read_cluster[mini_sign + 1][2]
        read_cluster[mini_sign] = [cluster_len, cluster_read_num, cluster_SV]
        del read_cluster[mini_sign + 1]

    read_cluster.sort(key=lambda x: x[1])
    if len(read_cluster)==2:
        yuzhi=min(abs(read_cluster[0][0]),abs(read_cluster[1][0]))*0.2
        if abs(read_cluster[0][0]-read_cluster[1][0])<=yuzhi:
            cluster_read_num = read_cluster[0][1] + read_cluster[1][1]
            cluster_len = (read_cluster[0][0] * read_cluster[0][1] + \
                           read_cluster[1][0] * read_cluster[1][1]) / cluster_read_num
            cluster_SV = read_cluster[0][2] + read_cluster[1][2]
            read_cluster[0] = [cluster_len, cluster_read_num, cluster_SV]
            del read_cluster[1]

    for cluster in read_cluster:
        average_len=cluster[0]+5
        #average_cov = sum(coverage) / len(coverage)
        if cluster[1]>=3 and cluster[1]>(0.2*read_num):
            #coverage.append(read_num)
            single_read_SV = cluster[2][0]
            for clutster_read in cluster[2]:
                if abs(clutster_read[0] - average_len)<max(0.1*average_len,10):
                    if calculate_SV_num(clutster_read[1]) < calculate_SV_num(single_read_SV[1]):
                        single_read_SV = clutster_read
                    elif calculate_SV_num(clutster_read[1]) == calculate_SV_num(single_read_SV[1]) and \
                            abs(clutster_read[0] - average_len) < abs(single_read_SV[0] - average_len):
                        single_read_SV = clutster_read
            if single_read_SV[1]!=[]:
                chr = single_read_SV[1][0][0]
                SV_pos1 = single_read_SV[1][0][1]
                SV_pos2 = single_read_SV[1][0][2]
                type = str(single_read_SV[1][0][4])

                if chr in main_chr:
                    ref = genes[chr][SV_pos1:SV_pos2]
                    if check_ref(ref) == True:
                        for single_sv in single_read_SV[1][1:]:
                            type = type + '+' + str(single_sv[4])
                            if single_sv[1] < SV_pos1:
                                SV_pos1 = single_sv[1]
                            if single_sv[2] > SV_pos2:
                                SV_pos2 = single_sv[2]
                        if abs(average_len) >= 40 and abs(average_len) < 50000:
                            FindCSV.append([chr,SV_pos1,SV_pos2,average_len,type,
                                    single_read_SV[1]])
                            #print(type)
                            '''w.write(
                                str(chr) + ' ' + str(SV_pos1) + ' ' + str(SV_pos2) + ' ' + str(
                                    average_len) + ' ' + type + ' ' + str(
                                    single_read_SV[1]) + '\n')'''
                        elif type == 'INV':
                            FindCSV.append([chr, SV_pos1, SV_pos2, single_read_SV[1][0][3], type,
                                single_read_SV[1]])
                            '''w.write(
                                str(chr) + ' ' + str(SV_pos1) + ' ' + str(SV_pos2) + ' ' + str(
                                    single_read_SV[1][0][3]) + ' ' + type + ' ' + str(
                                    single_read_SV[1]) + '\n')'''

###how to change the SV into figure
figure_save_path = "SV_image"
if not os.path.exists(figure_save_path):
    os.makedirs(figure_save_path)

for site in FindCSV:
    chr_name = site[0]
    pos1 = int(site[1])
    pos2 = int(site[2])
    SV_len = pos2-pos1+1
    # set the how many I replace a M
    num_i_to_m = int(SV_len / 20)
    num_i_to_m = max(10, num_i_to_m)
    # matrix length, the matrix is larger than SV, max flank length is 500
    matrix_len = SV_len + 2 * min(500, 2 * SV_len)
    Image_len = SV_len + 2 * min(500, 2 * SV_len)
    # Iamge name
    Image_name = str(site[0]) + '_' + str(site[1]) + '-' + str(site[2]) + ".png"
    # WIN_S is the start pos of window
    WIN_S = pos1 - min(500, 2 * SV_len)
    # if a SV is too long, the read number is too small, so we need to expand read
    read_expand_time = int(SV_len / 1000 + 1)
    # read_expand_time=100
    read_cigar = read_cigarfile(samfile, genes, chr_name, pos1, pos2)
    total_read = list(set(sorted(read_cigar.keys())))
    if len(total_read) != 0:
        # print(read_cigar)
        matrix = create_matrix(read_cigar, total_read)
        # print(matrix)
        create_image(matrix, Image_name, read_expand_time)
        # print(matrix)


localtime = time.strftime('%Y-%m-%d',time.localtime(time.time()))
w= open('FindCSV_result.vcf', 'w')
w.write('##fileformat=VCFv4.2'+'\n')
w.write('##source=FindCSV'+'\n')
w.write('##fileDate='+str(localtime)+'\n')

w.write('##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">'+'\n')
w.write('##ALT=<ID=DEL,Description="Deletion relative to the reference">'+'\n')
w.write('##ALT=<ID=INV,Description="Inversion of reference sequence">'+'\n')
w.write('##ALT=<ID=BND,Description="Breakend of translocation">'+'\n')

for chrom in genes.keys():
    w.write('##contig=<ID='+str(chrom) +',length='+str(len(genes[chrom]))+'>'+'\n')

w.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'+'\n')
w.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'+'\n')
w.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'+'\n')
w.write('##INFO=<ID=BP,Number=A,Type=String,Description="Breakpoint information of the variant described in this record">'+'\n')
w.write('##FILTER=<ID=q5,Description="Quality below 5">'+'\n')
w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'+'\n')

w.write('##CommandLine=FindCSV'+str(args.bam)+' '+str(args.fasta)+'\n')
w.write('#CHROM    POS ID REF ALT QUAL FILTER INFO FORMAT NULL'+'\n')



del_IDnum=0
ins_IDnum=0
inv_IDnum=0
csv_IDnum=0
for site in FindCSV:
    if str(site[4])=='DEL':
        del_IDnum+=1
        ID='FindCSV.del.'+str(del_IDnum)
    elif str(site[4])=='INS':
        ins_IDnum+=1
        ID = 'FindCSV.ins.' + str(ins_IDnum)
    elif str(site[4])=='INV':
        inv_IDnum+=1
        ID = 'FindCSV.inv.' + str(inv_IDnum)
    else:
        inv_IDnum += 1
        ID = 'FindCSV.inv.' + str(inv_IDnum)

    w.write(str(site[0]) + ' ' + str(site[1]) + ' '+ ID+' '+str(genes[str(site[0])][site[1]:site[2]+1])+' . 60 PASS ')
    w.write('SVTYPE='+str(site[4])+';'+'SVLEN='+str(int(site[3]))+';'+'END='+str(site[2])+';BP=')
    for subsite in site[5]:
        w.write(str(site[0])+':'+str(site[1])+'-'+str(site[2])+'|')
    w.write('; GT:DR:DV:PL:GQ   ./.:.:12:.,.,.:.'+'\n')
#I add modified tep to modify the DEL pos. and add soft clip read step.

localtime = time.asctime(time.localtime(time.time()))
print ('End',localtime)