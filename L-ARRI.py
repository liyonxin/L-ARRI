import os
import itertools

### Part1
os.system('mkdir arg-mge Nanofilt_result all_centrifuge raw_mges extract_mges mges_abundance 3_hbp_abundance ARRI')
os.system("mkdir ARG_raw_result ARG_abundance ARG_centrifuge ARG_reads extract_ARG arg_hbp")
path=os.getcwd()
sample_lst=[]
sample=os.listdir(path+'/rawdata')
lst_0=[]
outfile0=open(path+'/sum_length.txt','w')
for item in sample:
    line=item[0:item.index('.fastq.gz')] 
    sample_lst.append(line)    
for item in sample_lst:
    os.system('seqkit stats %s -T -o %s_length.txt'%(path+'/rawdata/%s.fastq.gz'%item,item))
lst_0=[]
for item in sample_lst:
    with open(path+'/%s_length.txt'%item,'r') as file1:
        f1=file1.readlines()
        line=f1[1].split('\t')
        lst_0.append(line[4].strip())
        outfile0.writelines('%s\t%s\n'%(item,line[4]))
    os.system('rm -r %s'%(path+'/%s_length.txt'%item))
outfile0.close()
### 0
for item0 in sample_lst:    
	os.system('gunzip -c %s | chopper -q 10 -l 500 | gzip > %s '%(path+'/rawdata/%s.fastq.gz'%item0, path+'/Nanofilt_result/%s_nanofilt.gz'%item0))
	os.system('seqtk seq -A  %s > %s'%(path+'/Nanofilt_result/%s_nanofilt.gz'%item0,path+'/Nanofilt_result/%s.fa'%item0) )

### 1
for item0 in sample_lst:
    os.system('centrifuge -f  -x /home/liyx/centrifuge/p+h+v/p_compressed+h+v  -U %s --report-file %s -S %s -p 56'%(path+'/Nanofilt_result/%s.fa'%item0,path+'/all_centrifuge/%s_all_report.tsv'%item0,path+'/all_centrifuge/%s_all_result.tsv'%item0))

### 2   
for item0 in sample_lst:
    os.system('minimap2 -x map-ont   -t  56 /home/liyx/minimap2/SARG_20211207_14210_filter.ffn  %s  > %s '%(path+'/Nanofilt_result/%s.fa'%item0,path+'/ARG_raw_result/%s_ARG'%item0))

### 3
os.system('lastdb -P56 -q -c trandb /home/liyx/database/mges/mobile-OG/mobileOG-db_beatrix-1.6.All.faa')
for item0 in sample_lst:    
    os.system('last-train -P56 --codon trandb %s > %s '%(path+'/Nanofilt_result/%s.fa'%item0,'%s_codon.train'%item0))
    os.system('lastal -P56 -p %s -m100 -D1e9 -K1 trandb %s > %s'%('%s_codon.train'%item0,path+'/Nanofilt_result/%s.fa'%item0,'%s_out.maf'%item0) )
    os.system('maf-convert psl %s > %s'%('%s_out.maf'%item0,path+'/raw_mges/%s_my-alignments.psl'%item0) )
    os.system('mv %s_codon.train %s'%(item0,path+'/raw_mges'))
    os.system('mv %s_out.maf %s'%(item0,path+'/raw_mges'))

for item0 in sample_lst:
    outfile1=open(path+'/extract_ARG/%s_ARG_cov_identity.txt'%item0,'w')
    with open(path+'/ARG_raw_result/%s_ARG'%item0,'r') as file1:
        f1=file1.readlines()
    lst1=[]
    for i in range(0,len(f1)):
        F1=f1[i].split('\t')
        lst1.append(F1)
    for i2 in range(0,len(lst1)):
        if len(lst1[i2])>=18:
            if float((float(lst1[i2][8])-float(lst1[i2][7]))/float(lst1[i2][6]))>0.9:
                if float(lst1[i2][9])/(float(lst1[i2][8])-float(lst1[i2][7]))>0.75:
                    outfile1.writelines(f1[i2])
    outfile1.close()

#### 4-2. ARGid-Abundance
for item0 in sample_lst:
    outfile2=open(path+'/ARG_abundance/%s_abundance_ARG.txt'%item0,'w')
    X=int(sample_lst.index('%s'%item0))
    k=float(lst_0[X])  
    list1=[]
    with open(path+'/extract_ARG/%s_ARG_cov_identity.txt'%item0,'r') as F:
        file1=F.readlines()
    for i in range(0,len(file1)):
        x=file1[i].split('\t')
        list1.append(x)
    list2=[]
    for item3 in list1:
        if item3[5] not in list2:
            list2.append(item3[5])
    if len(list2)>0:
        for item4 in list2:
            sum=0
            for i2 in range(len(list1)):        
                if list1[i2][5] ==item4:
                    x2=float((float(list1[i2][8])-float(list1[i2][7]))/float(list1[i2][6]))
                    sum=sum+x2
            ARGi_abundance=sum/(k/10**9)  
            outfile2.writelines("%s\t%s\n"%(item4,ARGi_abundance))
    outfile2.close()



#### 6. Extract-MGE
for item0 in sample_lst:
    outfile_2mges1=open(path+'/extract_mges/%s_same-reads_result.txt'%item0,'w')
    outfile_1mges1=open(path+'/extract_mges/%s_cov-idty.txt'%item0,'w')
    lst_mges1=[]  
    with open(path+'/raw_mges/%s_my-alignments.psl'%item0,'r') as f:
        file_1mges1 = f.readlines()
        for i in range(0,len(file_1mges1)):
            F1=file_1mges1[i].split('\t')
            lst_mges1.append(F1) 
    lst_2mges0=[]
    for i2 in range(0,len(lst_mges1)):
        if float((float(lst_mges1[i2][12])-float(lst_mges1[i2][11]))/float(lst_mges1[i2][10]))>0.7:
            if float(lst_mges1[i2][0])/(float(lst_mges1[i2][12])-float(lst_mges1[i2][11]))>0.7:
                outfile_1mges1.writelines(file_1mges1[i2])
                lst_2mges0.append(file_1mges1[i2]) 
    lst_2mges1=[]   
    for i in range(0,len(lst_2mges0)):
        F1=lst_2mges0[i].split('\t')
        lst_2mges1.append(F1)
    lst_2mges2=[]   
    for i in range(len(lst_2mges1)):
        item=lst_2mges1[i][13]
        if item not in lst_2mges2:
            lst_2mges2.append(item) 
##############################     
    for i in range(len(lst_2mges2)):
        lst_2mges3=[] 
        for i2 in range(len(lst_2mges1)):                
            if lst_2mges2[i]==lst_2mges1[i2][13]:
                lst_2mges3.append(lst_2mges1[i2])
        nrow=len(lst_2mges3)   
        lst_2mges4=[]
        combinations_lst = list(itertools.combinations(lst_2mges3, 2)) 
        for i3 in range(len(combinations_lst)):
            x=max(float(combinations_lst[i3][0][15]),float(combinations_lst[i3][1][15]))
            y=min(float(combinations_lst[i3][0][16]),float(combinations_lst[i3][1][16]))
            a=float(combinations_lst[i3][0][15])
            b=float(combinations_lst[i3][0][16])
            if (y-x)/(b-a)>0.8:
                if combinations_lst[i3][0] not in lst_2mges4:
                    lst_2mges4.append(combinations_lst[i3][0])
                if  combinations_lst[i3][1] not in lst_2mges4:
                    lst_2mges4.append(combinations_lst[i3][1])
        for i_2 in range(len(lst_2mges3)):
            if lst_2mges3[i_2] not in lst_2mges4:
                outfile_2mges1.writelines('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(lst_2mges3[i_2][0],lst_2mges3[i_2][8],lst_2mges3[i_2][9],lst_2mges3[i_2][10],lst_2mges3[i_2][11],lst_2mges3[i_2][12],lst_2mges3[i_2][13],lst_2mges3[i_2][14],lst_2mges3[i_2][15],lst_2mges3[i_2][16],lst_2mges3[i_2][17]))            
        sorted_lst = sorted(lst_2mges4, key=lambda x: float(x[15])) 
        lst_2mges5=[]
        if len(lst_2mges5)>0:
            for i4 in range(len(sorted_lst)-1):  
                x1=float(sorted_lst[i4][15])
                x2=float(sorted_lst[i4+1][15])
                length1=float(sorted_lst[i4][10])       
                if x2-x1<length1*2/5:
                    if sorted_lst[i4] not in lst_2mges5:
                        lst_2mges5.append(sorted_lst[i4])
                    if sorted_lst[i4+1] not in lst_2mges5:
                        lst_2mges5.append(sorted_lst[i4+1])
                else:
                    identity_lst=[]
                    for i5 in range(len(lst_2mges5)):
                        identity=float(lst_2mges5[i5][0])/float(lst_2mges5[i5][10])
                        identity_lst.append(identity)
                    identity_max_index=identity_lst.index(max(identity_lst))
                    result=lst_2mges5[identity_max_index]
                    for i7 in range(len(result)):
                        result[i7]=result[i7].strip()
                    outfile_2mges1.writelines('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(result[0],result[8],result[9],result[10],result[11],result[12],result[13],result[14],result[15],result[16],result[17]))
                    lst_2mges5=[]
            identity_lst2=[]
            for i6 in range(len(lst_2mges5)):
                identity=float(lst_2mges5[i6][0])/float(lst_2mges5[i6][10])
                identity_lst2.append(identity)
                identity_max_index=identity_lst2.index(max(identity_lst2))
                result2=lst_2mges5[identity_max_index]
            outfile_2mges1.writelines('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(result2[0],result2[8],result2[9],result2[10],result2[11],result2[12],result2[13],result2[14],result2[15],result2[16],result2[17]))
    outfile_2mges1.close()
    outfile_1mges1.close()


#### 7. MGE-Abundance
for item0 in sample_lst:
    outfile_3mges1=open(path+'/mges_abundance/%s_abundance_mges.txt'%item0,'w')
    list_3mges1=[]
    with open(path+'/extract_mges/%s_same-reads_result.txt'%item0,'r') as F:
        file_3mges1=F.readlines()
    for i in range(0,len(file_3mges1)):
        x=file_3mges1[i].split('\t')
        list_3mges1.append(x) 
    list_3mges2=[]     
    for item_3mges1 in list_3mges1:
        if item_3mges1[2] not in list_3mges2:
            list_3mges2.append(item_3mges1[2])
    for item_3mges2 in list_3mges2:
        sum=0 
        for i2 in range(len(list_3mges1)):        
            if list_3mges1[i2][2]==item_3mges2:
                x2=float((float(list_3mges1[i2][5])-float(list_3mges1[i2][4]))/float(list_3mges1[i2][3]))
                sum=sum+x2
        X=int(sample_lst.index('%s'%item0))
        k=float(lst_0[X])   
        MGE_abundance=sum/(k/10**9)  
        outfile_3mges1.writelines("%s\t%s\n"%(item_3mges2,MGE_abundance))     
    outfile_3mges1.close()

#### 8. ARG-MGE-Abundance
for item0 in sample_lst:
    lst0=[]
    outfile1=open(path+'/arg-mge/%s_arg-mge-result.txt'%item0,'w')
    outfile2=open(path+'/arg-mge/%s_arg-mge_abundance.txt'%item0,'w')
    file1=open(path+'/extract_ARG/%s_ARG_cov_identity.txt'%item0,'r')    
    lst_1=[]
    for line in file1.readlines():
        line=line.strip()
        lst_1.append(line)
    lst1=[]
    for item in lst_1:
        item=item.split('\t')
        lst1.append(item)
    lst_2=[]
    lst2=[]
    file2=open(path+'/extract_mges/%s_same-reads_result.txt'%item0,'r')
    for line in file2.readlines():
        line=line.strip()
        lst_2.append(line)
    for item in lst_2:
        item=item.split('\t')
        lst2.append(item)
    for i in range(len(lst1)):
        for i2 in range(len(lst2)):
            if lst1[i][0]==lst2[i2][6]:
                outfile1.writelines('%s\t%s\t%s\n'%(lst2[i2][2],lst1[i][5],lst2[i2][6]))
                lst0.append(lst2[i2][6])
    lst3=[] 
    for i in range(len(lst0)):
        if lst0[i] not in lst3:
            lst3.append(lst0[i])
    for item in lst3:  
        sum=0
        X=int(sample_lst.index('%s'%item0))
        k=float(lst_0[X])
        for i in range(len(lst1)):        
            if lst1[i][0]==item:
                sum=sum+float(lst1[i][9])/float(lst1[i][6]) 
        abundance=sum/(k/10**9)    
        outfile2.writelines('%s\t%s\n'%(item,abundance))
    outfile1.close()
    outfile2.close()

#### 9. HBP-Abundance  ARG-HBP-Abundance  ARG-MGE-HBP-Abundance  
for item0 in sample_lst: 
    outfile1=open(path+'/3_hbp_abundance/%s_HBP_abundance.txt'%item0,'w')
    outfile2=open(path+'/3_hbp_abundance/%s_ARG-HBP_abundance.txt'%item0,'w') 
    outfile3=open(path+'/3_hbp_abundance/%s_ARG-MGE-HBP_abundance.txt'%item0,'w')
    m=sample_lst.index(item0)
    k=float(lst_0[int(m)]) 
    with open(path+'/new-species-all.txt','r') as file1:
        f1=file1.readlines()
    lst1=[]
    for item1 in f1:
        line=item1.strip()
        line=line.split('\t')
        lst1.append(line) 
    with open(path+'/all_centrifuge/%s_all_result.tsv'%item0,'r') as file2:
        f2=file2.readlines()
    lst2=[]
    for item2 in f2:
        line=item2.strip()
        line=line.split('\t')
        if line[7]=='1':
            lst2.append(line)
    with open(path+'/all_centrifuge/%s_all_report.tsv'%item0,'r') as file3:
        f3=file3.readlines()
    lst3=[]
    for item3 in f3:
        line=item3.strip()
        line=line.split('\t')
        lst3.append(line)
    lst4=[]  
    for i in range(len(lst1)):
        sum=0
        for i2 in range(len(lst2)):
            if lst1[i][1]==lst2[i2][2]:
                lst4.append(lst2[i2])  
                sum=sum+float(lst2[i2][6])
        for i3 in range(len(lst3)):
            if lst1[i][1]==lst3[i3][1]:
                M=float(lst3[i3][3])
        if sum>0:
            if M>0:
                N=(sum/M)/(k/10**9)
                outfile1.writelines('%s\t%s\n'%(lst1[i][0],N))
    with open (path+'/extract_ARG/%s_ARG_cov_identity.txt'%item0,'r') as file4:
        f4=file4.readlines()
    lst5=[] 
    for i in range(len(f4)):
        line=f4[i].strip()
        line=line.split('\t')
        lst5.append(line)  
    for i in range(len(lst4)): 
        sum2=0
        for i2 in range(len(lst5)):
            if lst4[i][0]==lst5[i2][0]:
                sum2=sum2+(float(lst5[i2][8])-float(lst5[i2][7]))/float(lst5[i2][6])
        if sum2>0:
            ARG_HBP_abundance=sum2/(k/10**9) 
            outfile2.writelines('%s\t%s\n'%(lst4[i][0],ARG_HBP_abundance))
    lst0=[] 
    lst8=[] 
    lst6=[]
    file6=open(path+'/extract_mges/%s_same-reads_result.txt'%item0,'r')
    for line in file6.readlines():
        line=line.strip()
        line=line.split('\t')
        lst6.append(line)
    for i in range(len(lst5)):
        for i2 in range(len(lst6)):
            if lst5[i][0]==lst6[i2][6]:
                if lst6[i2][6] not in lst8:
                    lst8.append(lst6[i2][6])  
    for i in range(len(lst4)):
        for i2 in range(len(lst8)):
            if lst4[i][0]==lst8[i2]:
                lst0.append(lst8[i2])
    for item in lst0:
        sum3=0
        for i in range(len(lst5)):        
            if lst5[i][0]==item:
                sum3=sum3+(float(lst5[i2][8])-float(lst5[i2][7]))/float(lst5[i2][6])
        if sum3>0:
            ARG_MGE_HBP_abundance=sum3/(k/10**9) 
            outfile3.writelines('%s\t%s\n'%(item,ARG_MGE_HBP_abundance))
    outfile1.close()
    outfile2.close()
    outfile3.close()

#### 10. calculate-ARRI
for item0 in sample_lst:
    x=int(sample_lst.index('%s'%item0))
    k=float(lst_0[x])
    with open(path+'/ARG_abundance/%s_abundance_ARG.txt'%item0,'r') as file1:
        f1=file1.readlines()
    lst1=[]
    sum1=0
    for i in range(len(f1)):
        line=f1[i].split('\t')
        sum1=sum1+float(line[1].strip())
    lst1.append(item0)
    lst1.append(sum1)
    with open(path+'/mges_abundance/%s_abundance_mges.txt'%item0,'r') as file2:
        f2=file2.readlines()
    lst2=[]
    sum2=0
    for i2 in range(len(f2)):
        line=f2[i2].split('\t')
        sum2=sum2+float(line[1].strip())
    lst2.append(item0)
    lst2.append(sum2)
    with open(path+'/3_hbp_abundance/%s_HBP_abundance.txt'%item0,'r') as file3:
        f3=file3.readlines()
    lst3=[]
    lst_centrifuge=[]
    sum3=0
    for i3 in range(len(f3)):
        line=f3[i3].split('\t')
        sum3=sum3+float(line[1].strip())
    lst3.append(item0)
    lst3.append(sum3)
    with open(path+'/3_hbp_abundance/%s_ARG-HBP_abundance.txt'%item0,'r') as file4:
        f4=file4.readlines()
    lst4=[]
    sum4=0
    for i4 in range(len(f4)):
        line=f4[i4].split('\t')
        sum4=sum4+float(line[1].strip())
    lst4.append(item0)
    lst4.append(sum4)
    with open(path+'/3_hbp_abundance/%s_ARG-MGE-HBP_abundance.txt'%item0,'r') as file5:
        f5=file5.readlines()
    lst5=[]
    sum5=0
    for i5 in range(len(f5)):
        line=f5[i5].split('\t')
        sum5=sum5+float(line[1].strip())
    lst5.append(item0)
    lst5.append(sum5)             
    with open(path+'/arg-mge/%s_arg-mge_abundance.txt'%item0,'r') as file6:
        f6=file6.readlines()
    lst6=[]
    sum6=0
    for i6 in range(len(f6)):
        line=f6[i6].split('\t')
        sum6=sum6+float(line[1].strip())
    lst6.append(item0)
    lst6.append(sum6)
    outfile1=open(path+'/ARRI/%s_ARRI.txt'%item0,'w')
    ARRI=int((float(lst4[1])+float(lst5[1])+float(lst6[1]))/(float(lst1[1])+float(lst2[1])+float(lst3[1]))*10**4)
    outfile1.writelines('%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s\n%s\t%s'%('ARG-MGE',lst6[1],'ARG-HBP',lst4[1],'ARG-MGE-HBP',lst5[1],'Abundance(ARGs)',lst1[1],'Abundance(MGEs)',lst2[1],'Abundance(HBPs)',lst3[1],'ARRI',ARRI))
    outfile1.close()

#### 5. subARG/ARG-Abundance 
with open(path+'/structure_20181107.LIST','r') as file00:
    f00=file00.readlines()   
lst0=[]
for line in f00:
    line=line.strip()
    line=line.split('\t') 
    line[1]=line[1].replace('[','')
    line[1]=line[1].replace(']','')
    line[1]=line[1].replace("'",'')
    line[1]=line[1].replace(' ','')
    line[1]=line[1].split(',')
    lst0.append([line[0],line[1]])

for item0 in sample_lst:
    X1=int(sample_lst.index('%s'%item0))
    k=float(lst_0[X1]) 
    X=k/10**9
    outfile1=open(path+'/ARG_abundance/%s_subARG_abundance2.txt'%item0,'w')
    outfile2=open(path+'/ARG_abundance/%s_ARG_abundance3.txt'%item0,'w')
    LST1=[]    
    with open(path+'/extract_ARG/%s_ARG_cov_identity.txt'%item0,'r') as file1:
        f1=file1.readlines()    
    for line in f1:
        line=line.split('\t')
        line[13]=line[13].split(':')
        for item1 in lst0:
            if line[13][0] in item1[1]:
                LST1.append([line[0],item1[0],(float(line[16])-float(line[15]))/float(line[14])/X])
    LST2=[]  
    lst2=[]
    for item in LST1:
        if item[1] not in lst2:
            lst2.append(item[1])
    for item1 in lst2:
        sum=0
        for item2 in LST1:
            if item2[1]==item1:
                sum=sum+float(item2[2])
        LST2.append([item1,sum])
        outfile1.writelines('%s\t%s\n'%(item1,sum))
    lst3=[]
    for item in LST2:
        line=item[0].split('_')
        line=line[0] 
        if line not in lst3:
            lst3.append(line)
    for item1 in lst3:
        sum2=0
        for item2 in LST2:
            if item1 in item2[0] :
                sum2=sum2+float(item2[1])
        outfile2.writelines('%s\t%s\n'%(item1,sum2))
    outfile1.close()
    outfile2.close()
