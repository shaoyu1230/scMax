reportdir=$1
#reportdir需要是CLUSTER和ANNO的输出目录，report.html也会生成在同级目录
cp /annogene/data2/bioinfo/SCV/PROJECT/RD_project/SCV_RD_2503/bin/report.ipynb ${reportdir}
cd ${reportdir}
/annogene/data2/bioinfo/PMO/wenzhudai/Software/mamba/envs/jupyter_dwz/bin/jupyter nbconvert --no-input --template pj --to html report.ipynb
zip -r report.zip 1_cluster_characterization  2_marker_expression  3_annotation  4_celltype_fraction  5_celltype_characterization report.html
