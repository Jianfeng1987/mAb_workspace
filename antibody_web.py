# app.py
import streamlit as st
import subprocess
import os
import pandas as pd
from pathlib import Path

# ======================
# 页面配置
# ======================
st.set_page_config(
    page_title="抗体可变区分析 & 进化树",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 抗体可变区序列提取 + 进化树分析")
st.markdown("上传 FASTA/FASTQ 文件，自动完成 IgBLAST → CDR 提取 → 配对 → 进化树构建。")

# ======================
# 侧边栏参数配置
# ======================
st.sidebar.header("⚙️ 分析参数")
species = st.sidebar.selectbox(
    "物种",
    ["mouse", "human", "rabbit", "rat"],
    index=0
)
cdr_def = st.sidebar.selectbox(
    "CDR 定义",
    ["kabat", "imgt", "chothia", "contact", "martin"],
    index=0
)
run_pair = st.sidebar.checkbox("进行 VH/VL 配对", value=True)
run_tree = st.sidebar.checkbox("构建进化树", value=True)

# ======================
# 【新增】侧边栏使用说明（含命名格式）
# ======================
st.sidebar.markdown("---")
st.sidebar.markdown("### 📖 使用说明")
st.sidebar.markdown(
    """
    1. **准备文件**：确保序列文件命名规范。
    2. 上传测序文件（FASTA / FASTQ）。
    3. 选择物种与 CDR 定义。
    4. 勾选是否配对 / 进化树。
    5. 点击「开始分析」。
    6. 下载 Excel 与进化树文件。
    """
)

st.sidebar.markdown("### 📝 FASTA 命名格式")
st.sidebar.markdown(
    """
    **必须包含孔位信息**，否则无法配对！
    
    **格式：**
    `<样本>-B<板号>-<孔位>-[H/L]-<序号>`
    
    **示例：**
    - `80-2-B1-96A01-H-001` (重链)
    - `80-2-B1-96A01-L-001` (轻链)
    
    **说明：**
    - `-H-` 代表重链，`-L-` 代表轻链。
    - `96A01` 对应 96 孔板的 A1 孔。
    """
)

# ======================
# 文件上传
# ======================
uploaded_file = st.file_uploader(
    "📄 上传序列文件 (FASTA / FASTQ)",
    type=["fasta", "fa", "fastq", "fq"]
)

# 【新增】文件上传后的即时提示
if uploaded_file is not None:
    st.info(f"✅ 已接收文件：{uploaded_file.name}。请确保文件名包含孔位信息（如 96A01），否则配对功能将失效。")
else:
    st.warning("⚠️ 请上传 FASTA/FASTQ 文件。注意：文件名需包含孔位编号（如 96A01）以支持 VH/VL 配对。")

# ======================
# 主逻辑
# ======================
if uploaded_file is not None:
    # 保存上传文件
    input_path = Path("input.fasta")
    input_path.write_bytes(uploaded_file.read())

    if st.button("🚀 开始分析", type="primary"):
        with st.spinner("分析中，请稍候…（IgBLAST 和 MUSCLE 可能较慢）"):
            # 构造命令
            cmd = [
                "python3",
                "antibody_vr_extractor_CDR_final.py",
                "--input", str(input_path),
                "--species", species,
                "--cdr_def", cdr_def,
                "--skip_setup"
            ]
            if run_pair:
                cmd.append("--pair")
            if run_tree:
                cmd.append("--tree")

            # 执行脚本
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=os.getcwd()
            )

            if result.returncode != 0:
                st.error("❌ 分析失败")
                st.code(result.stderr, language="bash")
            else:
                st.success("🎉 分析完成！")

                # ======================
                # 结果展示与下载
                # ======================
                work_dir = Path("antibody_work")

                # 1. Excel 结果
                excel_files = list(work_dir.glob("*_VH_VL_paired.xlsx"))
                if not excel_files:
                    excel_files = list(work_dir.glob("*.xlsx"))

                if excel_files:
                    excel_path = excel_files[0]
                    st.markdown("### 📊 分析结果 Excel")
                    with open(excel_path, "rb") as f:
                        st.download_button(
                            label="⬇️ 下载 Excel 结果",
                            data=f,
                            file_name=excel_path.name,
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                        )
                    
                    # 预览前 5 行
                    try:
                        df = pd.read_excel(excel_path, sheet_name=0)
                        st.dataframe(df.head())
                    except Exception:
                        pass

                # 2. 进化树文件
                tree_dir = work_dir / "phylogenetic_trees"
                if tree_dir.exists():
                    st.markdown("### 🌳 进化树文件")

                    for ext in ["svg", "png", "nwk"]:
                        files = list(tree_dir.glob(f"*.{ext}"))
                        if files:
                            st.markdown(f"**{ext.upper()} 文件**")
                            for f in files:
                                with open(f, "rb") as fp:
                                    st.download_button(
                                        label=f"⬇️ {f.name}",
                                        data=fp,
                                        file_name=f.name,
                                        mime=(
                                            "image/svg+xml" if ext == "svg"
                                            else "image/png" if ext == "png"
                                            else "text/plain"
                                        )
                                    )

                    # 显示 SVG 树
                    svg_files = list(tree_dir.glob("*.svg"))
                    if svg_files:
                        st.image(str(svg_files[0]), caption="VH 进化树 (SVG)")

                # 3. 日志输出（可选）
                with st.expander("📜 查看运行日志"):
                    st.code(result.stdout, language="bash")

else:
    st.info("👆 请先上传 FASTA 或 FASTQ 文件")