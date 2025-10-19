import os
import sys
import pandas as pd
import json

jsonfile = sys.argv[1]
outfile = sys.argv[2]

def read_json(infile):
    with open(infile, 'r', encoding='utf-8') as file:
        data = json.load(file)
    '''
    antismash results are in  data['records']
    '''
    tmplist = []
    for res in data['records']:
        contig_id = res['id']
        contig_name = res['name']
        contig_description = res['description']
        
        region_cn = 1
        previous_start = None
        previous_end = None
        
        for region in res['modules']['antismash.detection.hmm_detection']['rule_results']['cds_by_protocluster']:
            '''
            type of product
            '''
            print(region[0]['qualifiers']['category'])
            category = ",".join(region[0]['qualifiers']['category'])
            product = ",".join(region[0]['qualifiers']['product'])
            
            '''
            whether in previous region
            '''
            
            region_loc = region[0]['location']
            coord_part = region_loc.split('[')[1].split(']')[0]
            start, end = coord_part.split(':')
            # 一定要做转换，否则可以比较不报错但是结果不正确
            start = int(start)
            end = int(end)
            if previous_start == None:
                regionid = 'r'+str(region_cn)
                previous_start,previous_end = start,end
            else:
                if start <= previous_end:
                    regionid = 'r'+str(region_cn)
                    if end < previous_end:
                        pass
                    else:
                        previous_end = end
                else:
                    # print(f"{start} > {previous_end}")
                    region_cn += 1
                    regionid = 'r'+str(region_cn)
                    previous_start,previous_end = start,end
                    
                
            strand = region_loc.split('(')[1].split(')')[0]
            tmplist.append((contig_id,contig_name,contig_description,regionid,start,end,strand,category,product))
    
    df = pd.DataFrame(tmplist)
    df.columns = ['contig_id','contig_name','contig_description','regionid','start','end','strand','category','product']
    
    return df

if __name__ == '__main__':
    outdf = read_json(jsonfile)
    outdf.to_csv(outfile,sep="\t",index=None)
