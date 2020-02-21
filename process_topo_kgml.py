# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 11:09:30 2018

@author: libingrui
"""
import pandas as pd
import os
def process_topo_kgml():
    outfile=open(r'D:\pathway_akk.txt','w')
    with open(r'D:\g_id2symb.txt','r') as infile:
        id_symb={}
        for lines in infile:
            lines=lines.strip()
            lines=lines.split('\t')
            id_symb[lines[0]]=lines[1]
    all_filename=os.listdir('D:\Topo_kgml')
    all_filename_list=[]
    for filename in all_filename:
        filename=filename.replace(' ','_')
        all_filename_list.append(filename)
    with open (r'D:\module.txt','r') as module_file:
        num=0
        module_dict={}
        for module_line in module_file:
            num+=1
            module_line.strip()
            module_line=module_line.replace('\t','')
            module_line=str(module_line).replace('\n','').split(',')
            for module_unit in module_line:
                module_dict[module_unit]=num            
    all_symb_value=[]
    all_symb_symb=[]        
    for filename in all_filename:
        unit_symb_value=[]
        unit_symb_symb=[]
        readpath='D:\Topo_kgml\\'+str(filename)
        with open(readpath,'r') as kgmlfile:
            symb_value={}
            n=0
            for lines in kgmlfile:
                lines=lines.replace('\n','').strip()
                lines=lines.split(' ')
                if len(lines)>2 and lines[2]!='NA"':
                    lines[2]=lines[2][:-1]
                    lines[1]=lines[1][1:]
                    if lines[1] in id_symb.keys():
                        if id_symb[lines[1]] in module_dict.keys():
                            name='mod'+str(module_dict[id_symb[lines[1]]])
                            if name not in symb_value.keys():
                                symb_value[name]=lines[2]
                            else:
                                symb_value[name]=float(symb_value[name])+float(lines[2])
                        else:
                            symb_value[id_symb[lines[1]]]=lines[2]
        for v in symb_value.values():
            unit_symb_value.append(v)
        for s in symb_value.keys():
            unit_symb_symb.append(s)
        all_symb_value.append(unit_symb_value)
        all_symb_symb.append(unit_symb_symb)
    symb_value_df=pd.DataFrame(all_symb_value,index=all_filename_list)
    symb_value_df=symb_value_df.transpose()
    symb_value_df.to_csv('D:\pathway_akk.csv')
    symb_symb_df=pd.DataFrame(all_symb_symb,index=all_filename_list)
    symb_symb_df=symb_symb_df.transpose()
    symb_symb_df.to_csv('D:\pathway_contents.csv')
    kgmlfile.close()
	module_file.close()
    infile.close()



process_topo_kgml()