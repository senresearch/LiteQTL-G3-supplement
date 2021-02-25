#!/bin/bash 

lscpu
nvidia-smi 
#top -u xiaoqihu -b -n 10 -d 1 & P1=$!
python time_tensorqtl.py & PIDPYTHON=$!
top -p $PIDPYTHON -b -d 1 & PIDTOP=$! 
wait $PIDPYTHON $PIDTOP


