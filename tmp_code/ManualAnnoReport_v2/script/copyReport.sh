reportdir=$1
outdir=$2
# 检查条件并拷贝文件
if [ -d "${outdir}/2_marker_expression" ] && [ -f "${outdir}/1_cluster_characterization/1_Group.UMAP.png" ]; then
    cp ${reportdir}/report.ipynb  ${outdir}/report.ipynb
elif [ ! -d "${outdir}/2_marker_expression" ] && [ -f "${outdir}/1_cluster_characterization/1_Group.UMAP.png" ]; then
    cp ${reportdir}/report_noref.ipynb  ${outdir}/report.ipynb
elif [ -d "${outdir}/2_marker_expression" ] && [ ! -f "${outdir}/1_cluster_characterization/1_Group.UMAP.png" ]; then
    cp ${reportdir}/report_nogroup.ipynb  ${outdir}/report.ipynb
elif [ ! -d "${outdir}/2_marker_expression" ] && [ ! -f "${outdir}/1_cluster_characterization/1_Group.UMAP.png" ]; then
    cp ${reportdir}/report_noref.nogroup.ipynb  ${outdir}/report.ipynb
fi
