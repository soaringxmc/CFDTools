输入
run_average.txt 控制文件
y0_slice.3400000.dat tecplot转换后的文本文件point格式，平均场
输出
test.dat（中间结果） 平均场，point格式，和y0_slice.3400000.dat平均场数据相同，只是格式不同
.avg.dat（最终结果）


_cycle.mcr 
输入
不同时刻的表面二进制转换成文本。变量包括
"X" "Y" "Z" "RHO" "P" "T" "TAU_WALL" "YPLUS" "YPLUS_WM" "XPLUS" "DIL" "P_AVG" "P_RMS" "U-X" "U-Y" "U-Z"
若x和y相同，则平均
输出
airfoil.dat


run_average.exe
功能
表面数据时间平均和展向平均


