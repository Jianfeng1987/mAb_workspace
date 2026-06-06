抗体可变区序列分析与进化树构建工具 
项目简介 
本工具是一套完整的抗体序列分析流水线，支持从原始FASTA/FASTQ文件到最终结果的全自动化处理，包含IgBLAST比对、CDR提取、VH/VL配对、进化树构建四大核心功能。同时提供命令行界面和Web可视化界面，满足不同场景的使用需求。
功能特性 
✅ 支持多物种分析：mouse/human/rabbit/rat
✅ 多CDR定义标准：kabat/imgt/chothia/contact/martin
✅ 自动完成IgBLAST序列比对与可变区提取
✅ VH/VL孔位配对分析（依赖规范命名）
✅ 进化树构建（NJ法+p-distance，支持SVG/PNG/Newick格式导出）
✅ 结果自动生成Excel报表（含质控状态、基因信息、CDR序列）
✅ 双模式运行：命令行批量处理 / Web界面交互操作
✅ 支持腾讯云COS存储（可选）
环境要求 
操作系统 
•	Linux（推荐Ubuntu 22.04+）
•	Windows（需通过WSL 2运行）
•	macOS（需安装Homebrew依赖）
基础依赖 
类型	版本要求	说明
Python	3.8+	主运行环境
系统工具	hmmer、muscle	序列比对依赖（需通过apt/brew安装）
IgBLAST	1.22.0	抗体基因比对引擎（需本地tar包）
物种数据库	鼠源VDJ数据库	需本地tar包（如mouse_gl_VDJ.tar）
快速安装 
1. 克隆仓库 
bash
bash
git clone https://github.com/your-repo/antibody-analysis.git
cd antibody-analysis
2. 安装系统依赖 
bash
bash
# Ubuntu/Debian
sudo apt update && sudo apt install -y hmmer muscle

# macOS（需先安装Homebrew）
brew install hmmer muscle
3. 创建虚拟环境（推荐） 
bash
bash
python3 -m venv venv
source venv/bin/activate  # Linux/macOS
# Windows: venv\Scripts\activate
4. 安装Python依赖 
bash
bash
pip install -r requirements.txt
若遇系统保护报错（externally-managed-environment），可加参数强制安装：pip install -r requirements.txt --break-system-packages
5. 准备本地文件 
将以下文件放入项目根目录：
•	ncbi-igblast-1.22.0-x64-linux.tar.gz（IgBLAST安装包）
•	mouse_gl_VDJ.tar（鼠源数据库，其他物种需自行准备）
使用说明 
模式一：命令行批量处理 
基础用法 
bash
bash
# 交互式运行（默认）
python antibody_vr_extractor_CDR_final.py /path/to/your_sequences.fasta

# 非交互式（指定参数）
python antibody_vr_extractor_CDR_final.py \
  --input /path/to/your_sequences.fasta \
  --species mouse \
  --cdr_def kabat \
  --skip_setup \
  --pair \
  --tree
常用参数 
参数	说明	默认值
--input	输入FASTA/FASTQ文件路径	必填
--species	物种选择	mouse
--cdr_def	CDR定义标准	kabat
--skip_setup	跳过环境安装（已配置时）	否
--pair	执行VH/VL配对	否
--tree	构建进化树	否
________________________________________
模式二：Web可视化界面（推荐） 
启动服务 
bash
bash
streamlit run app.py --server.port 8501
浏览器访问：http://localhost:8501
操作流程 
1.	上传文件：拖拽FASTA/FASTQ文件到上传区（需符合命名规范）
2.	参数配置：在侧边栏选择物种、CDR定义，勾选是否配对/进化树
3.	开始分析：点击「🚀 开始分析」按钮，等待进度完成
4.	下载结果：直接点击Excel、进化树文件的下载按钮
关键规范：FASTA文件命名 
⚠️ 配对功能依赖严格的序列ID命名，否则无法识别孔位！
命名格式 
纯文本
纯文本
<样本前缀>-B<板号>-<孔位编号>-[H/L]-<序号>
•	H：重链（VH）；L：轻链（VL）
•	孔位编号：对应96孔板坐标（如96A01=A1孔）
正确示例 
fasta
fasta
>80-2-B1-96A01-H-001
ATGGGATCAAGCTGACCCAGTCTCCA...
>80-2-B1-96A01-L-001
ATGGGATCAAGCTGACCCAGTCTCCA...
错误示例 
fasta
fasta
>seq1  # 无孔位信息，无法配对
>80-2-H-001  # 缺少板号和孔位，无法配对
输出结果说明 
目录结构 
纯文本
纯文本
antibody_work/
├── antibody_results_kabat.xlsx       # 主结果Excel（含PASS/FAIL sheet）
├── antibody_results_kabat_VH_VL_paired.xlsx  # 配对结果Excel
└── phylogenetic_trees/               # 进化树文件目录
    ├── VH_tree.nwk                   # VH链Newick格式
    ├── VH_tree.svg                   # VH链矢量图（推荐用于论文）
    ├── VH_tree.png                   # VH链位图
    └── VL_tree.*                     # VL链对应文件
Excel文件说明 
Sheet名称	内容
PASS_kabat	质控通过的序列（含CDR、基因信息）
FAIL_kabat	质控失败的序列（含失败原因）
配对总览	孔位配对状态统计（PAIRED_PASS/PAIRED_PARTIAL等）
PAIRED_PASS序列	成功配对的VH/VL完整信息
常见问题解决 
1. IgBLAST安装失败 
•	检查ncbi-igblast-1.22.0-x64-linux.tar.gz是否在项目根目录
•	确认系统架构为x64（32位系统不支持）
2. 系统保护报错（externally-managed-environment） 
•	优先使用虚拟环境（推荐）
•	临时解决：pip install ... --break-system-packages
3. 进化树构建失败 
•	检查muscle是否安装：muscle -version
•	确认PAIRED_PASS序列数≥3（少于3条无法构建树）
4. 配对结果为空 
•	检查FASTA序列ID是否符合命名规范（重点：包含孔位编号和H/L标识）
扩展功能 
腾讯云COS存储（可选） 
运行时选择「使用COS存储」，输入SecretId、SecretKey等信息，结果会自动上传至指定存储桶。
自定义进化树样式 
修改antibody_vr_extractor_CDR_final.py中的_draw_beautiful_tree函数，可调整颜色、标记、标签格式。
贡献指南 
欢迎提交Issue和PR！请确保：
•	代码符合PEP8规范
•	新增功能包含测试用例
•	更新对应文档说明
许可证 
本项目采用MIT许可证，详情见LICENSE文件。
________________________________________
技术支持：如有问题，请在GitHub仓库提交Issue，或联系项目维护者。

