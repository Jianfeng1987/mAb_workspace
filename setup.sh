#!/bin/bash
# setup.sh - 云端部署专用安装脚本

set -e  # 遇到错误则停止

echo "开始安装云端依赖..."

WORK_DIR="/workspace/antibody_work"
IGBLAST_DIR="/workspace/ncbi-igblast-1.22.0"
DB_DIR="${WORK_DIR}/igblast_db"

# 1. 创建必要的目录
mkdir -p ${WORK_DIR} ${DB_DIR}

# 2. 下载并安装 IgBLAST (从 NCBI 官方镜像)
echo "正在从 NCBI 下载 IgBLAST 1.22.0..."
if [ ! -f "${IGBLAST_DIR}/bin/igblastn" ]; then
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/1.22.0/ncbi-igblast-1.22.0-x64-linux.tar.gz
    tar -xzf ncbi-igblast-1.22.0-x64-linux.tar.gz -C /workspace/
    rm ncbi-igblast-1.22.0-x64-linux.tar.gz
    chmod +x ${IGBLAST_DIR}/bin/igblastn
    echo "  ✅ IgBLAST 安装完成。"
else
    echo "  ✅ IgBLAST 已存在。"
fi

# 3. 下载并安装鼠源 IMGT 数据库
echo "正在从 IMGT 下载鼠源数据库..."
if [ ! -f "${DB_DIR}/mouse_gl_V.nhr" ]; then
    cd ${DB_DIR}
    # 下载所有必要的种系序列文件
    imgt_base="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/"
    files=("IGHV.fasta" "IGHD.fasta" "IGHJ.fasta" "IGKV.fasta" "IGKJ.fasta" "IGLV.fasta" "IGLJ.fasta")
    for file in "${files[@]}"; do
        wget -q "${imgt_base}${file}"
    done
    echo "  ✅ 数据库原始文件下载完成。"
    # 简化处理：合并并创建BLAST数据库
    cat IGHV.fasta IGKV.fasta IGLV.fasta | sed '/^>/! s/[ .]//g' > mouse_gl_V
    cat IGHD.fasta | sed '/^>/! s/[ .]//g' > mouse_gl_D
    cat IGHJ.fasta IGKJ.fasta IGLJ.fasta | sed '/^>/! s/[ .]//g' > mouse_gl_J
    # 使用刚安装的 IgBLAST 工具建库
    ${IGBLAST_DIR}/bin/makeblastdb -in mouse_gl_V -dbtype nucl -parse_seqids -out mouse_gl_V
    ${IGBLAST_DIR}/bin/makeblastdb -in mouse_gl_D -dbtype nucl -parse_seqids -out mouse_gl_D
    ${IGBLAST_DIR}/bin/makeblastdb -in mouse_gl_J -dbtype nucl -parse_seqids -out mouse_gl_J
    echo "  ✅ 鼠源 BLAST 数据库构建完成。"
else
    echo "  ✅ 数据库已存在。"
fi

echo "所有云端依赖安装完毕！"
