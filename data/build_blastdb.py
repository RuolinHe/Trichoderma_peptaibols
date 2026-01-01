import pandas as pd
import os

# Fungi_flag=1
Fungi_flag=0 #Trichoderma
if Fungi_flag==1:
    # 读取 Excel 文件中的数据
    excel_file = '/storage/disk1/HRL/Project/SM_in_genome/Fungi/Fungi_genomes_clean.xlsx'
    sheet_name = 'high_quality'
    data = pd.read_excel(excel_file, sheet_name=sheet_name)

    # 获取第一列的字符串列表
    genome_ids = data.iloc[:, 0].tolist()

    # 循环处理每个字符串
    for genome_id in genome_ids:
        # 构建输入和输出路径
        input_path = f'/storage/disk1/HRL/Project/SM_in_genome/Fungi/ncbi_dataset/data/{genome_id}/{genome_id}.faa'
        output_path = f'/storage/disk2/HRL/Project/4-Peptaibols_NRPS/data/blastdb/Fungi/{genome_id}'
        
        # 如果输入路径存在，执行命令
        if os.path.exists(input_path):
            command = f'makeblastdb -in {input_path} -dbtype prot -out {output_path}'
            os.system(command)
        else:
            print(f"输入路径 {input_path} 不存在，跳过 {genome_id} 的处理")

else: #Trichoderma
    # 指定.fa文件所在的目录
    faa_directory = '/storage/disk2/HRL/Project/4-Peptaibols_NRPS/data/genome/galaxy/'

    # 遍历目录中所有的.fa文件
    for file_name in os.listdir(faa_directory):
        if file_name.endswith('.faa'):
            # 提取基因组ID
            genome_id = os.path.splitext(file_name)[0]

            # 构建输入和输出路径
            input_path = os.path.join(faa_directory, file_name)
            output_path = f'/storage/disk2/HRL/Project/4-Peptaibols_NRPS/data/blastdb/Trichoderma/{genome_id}'
            
            # 执行命令
            command = f'makeblastdb -in {input_path} -dbtype prot -out {output_path}'
            os.system(command)
