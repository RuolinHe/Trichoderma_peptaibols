import os

# 定义目录路径
taxonomy='Fungi'
# taxonomy='Trichoderma'
pdb_directory = f'/storage/disk2/HRL/Project/4-Peptaibols_NRPS/data/blastdb/{taxonomy}/'

# 遍历目录中所有的.pdb文件
for file_name in os.listdir(pdb_directory):
    if file_name.endswith('.pdb'):
        # 提取基因组ID
        genome_id = os.path.splitext(file_name)[0]
        
        # 构建输入和输出路径
        input_path = f'/storage/disk2/HRL/Project/4-Peptaibols_NRPS/data/Tqa.fasta'
        output_path = f'/storage/disk2/HRL/Project/4-Peptaibols_NRPS/data/blastp/{taxonomy}/{genome_id}.txt'
        db_path = os.path.join(pdb_directory, genome_id)
        
        # 执行命令
        command = f'blastp -query {input_path} -db {db_path} -outfmt "6 qseqid qacc sseqid sallseqid sgi sallgi sacc sallacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp scovhsp" -out {output_path} -num_threads 120'
        os.system(command)
