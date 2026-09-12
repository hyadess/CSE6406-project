# export PATH=/home/mis/SimPhy/bin:$PATH

python run_project_parallel.py \
  --output work/paper_subset_50x10_L1600 \
  --taxa 21 --loci 50 --replicates 10 \
  --lengths 1600 \
  --models gtr_g4 gtr hky_g4 hky jc \
  --seed 6406001 \
  --skip-ils-gate \
  --workers 12 --threads AUTO

# python run_project_parallel.py \
#   --output work/paper_subset_50x10_L200 \
#   --taxa 21 --loci 50 --replicates 10 \
#   --lengths 200 \
#   --models gtr_g4 gtr hky_g4 hky jc \
#   --seed 6406001 --skip-ils-gate \
#   --workers 12 --threads AUTO

# python run_project_parallel.py \
#   --output work/paper_subset_50x10_L800 \
#   --taxa 21 \
#   --loci 50 \
#   --replicates 10 \
#   --lengths 800 \
#   --models gtr_g4 gtr hky_g4 hky jc \
#   --seed 6406001 \
#   --skip-ils-gate \
#   --workers 8 \
#   --threads AUTO 
