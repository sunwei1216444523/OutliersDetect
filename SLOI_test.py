#-*- coding: utf-8 -*-
# 这里使用SLOM中的初始网格数据来测试本文算法的有效性

from __future__ import division
import numpy
import math

raw_data= [[-16,-9,-16,4,8,25,-2,20,-11,9],
    [-3,-1,9,-12,1,1,-1,-2,-4,-2],
    [14,1,11,2,-13,15,4,3,11,19],
    [-10,16,-11,-2,-10,-11,-17,4,8,-15],
    [-5,20,-11,4,-5,8,6,6,-2,-1],
    [15,10,-9,7,12,-9,-18,16,8,-6],
    [0,0,0,0,-21,-5,12,-15,-5,11],
    [0,0,0,0,5,6,1,1,-9,3],
    [0,8,0,0,-9,-8,-1,-2,9,5],
    [0,0,0,0,19,-1,-2,-7,-3,-12]]
"""
定义一个使用最大值-最小值归一化给定二维数组的函数
"""
def normFunction(paraArray):
    normArray = [[[]for ii in range(len(paraArray[0]))] for i in range(len(paraArray))]
    maxValueRow = []
    minValueRow = []
    for item in paraArray:
        maxValueRow.append(max(item))
        minValueRow.append(min(item))
    maxValue = max(maxValueRow)
    minValue = min(minValueRow)
    for i in range(len(paraArray)):
        for j in range(len(paraArray[0])):
            normArray[i][j] = round((paraArray[i][j]-minValue)/(maxValue-minValue)*(1-0.1)+0.1,3)
    return normArray

"""
计算目标对象的邻域距离
target 目标对象属性值
neighbours 目标对象周围邻居的对应坐标数组
"""
def getD(target,neighbours):
    Local_data_diff_array = []
    # 按照建立的邻居索引把对应的归一化后的值添加到集合中便于后面的计算
    for item in neighbours:
        Local_data_diff_array.append(abs(target-original_data[item[0]][item[1]]))
    result = (sum(Local_data_diff_array)-max(Local_data_diff_array))/(len(Local_data_diff_array)-1)
    return result

"""
计算每一个对象的信息熵
neighbours 目标对象周围邻居的对应坐标数组
"""
def getStable(neighbours):
    neighbourValue = []
    for item in neighbours:
        neighbourValue.append(original_data[item[0]][item[1]])
    sumValue = sum(neighbourValue)
    attributEntropy = 0
    for i in neighbourValue:
        i = i+0.00000000001
        rate  = i/(sumValue+0.00000000001*len(neighbourValue))
        attributEntropy += rate*math.log(1/rate,2)
    REntropy = attributEntropy/math.log(len(neighbours),2)
    return round(REntropy,3)
   
original_data = []
result = [[[]for x in range(10)] for i in range(10)]         # 存储SLOI值
result_stable = [[[]for x in range(10)] for i in range(10)]  # 存储局部区域稳定值
result_d = [[[]for x in range(10)] for i in range(10)]       # 存储离群度绝对值d

# 数据归一化
max_v = 0
min_v = 0
max_row = []
min_row = []
for item in raw_data:
    max_row.append(max(item))
    min_row.append(min(item))
max_v = max(max_row)
min_v = min(min_row)
for x in range(10):
    row_ = []
    for y in range(10):
        row_.append((raw_data[x][y]-min_v)/(max_v-min_v))
    original_data.append(row_)




# 循环遍历所有位置获取每个对象邻居的坐标,并存储在neighborTable中
neighbourTable = [[[]for ii in range(10)] for i in range(10)]
for x in range(10):
    for y in range(10):
        neighbour = []
        if(x == 0):
            if(y == 0):
                neighbour = [(0,1),(1,0),(1,1)]
            else:
                if(y == 9):
                    neighbour = [(0,8),(1,8),(1,9)]
                else:
                    neighbour = [(x,y-1),(x,y+1),(x+1,y),(x+1,y-1),(x+1,y+1)]
        else:
            if(x == 9):
                if(y == 0):
                    neighbour = [(9,1),(8,0),(8,1)]
                else:
                    if(y == 9):
                        neighbour = [(9,8),(8,8),(8,9)]
                    else:
                        neighbour = [(x,y-1),(x,y+1),(x-1,y),(x-1,y-1),(x-1,y+1)]
            else:
                if(y == 0):
                    neighbour = [(x-1,y),(x+1,y),(x-1,y+1),(x,y+1),(x+1,y+1)]
                else:
                    if(y==9):
                        neighbour = [(x-1,y),(x+1,y),(x-1,y-1),(x,y-1),(x+1,y-1)]
                    else:
                        neighbour = [(x-1,y-1),(x-1,y),(x-1,y+1),(x,y-1),(x,y+1),(x+1,y-1),(x+1,y),(x+1,y+1)]
        neighbourTable[x][y] = neighbour


#计算每个对象的邻域距离d
for x in range(10):
    for y in range(10):
        d = getD(original_data[x][y],neighbourTable[x][y])
        result_d[x][y] = round(d,3)


# 计算稳定性系数 ;
for x in range(10):
    for y in range(10):
        result_stable[x][y] = getStable(neighbourTable[x][y])

"""
把稳定性系数归一化
"""
norm_stable = normFunction(result_stable)
print("==============norm_stable=============")
for i in norm_stable:
    print(i)


"""
重新计算SLOI
"""
for i in range(10):
    for j in range(10):
        # 不再对稳定度进行再一次的归一化，直接使用norm_stable,且使用arcsin
        result[i][j] =round(result_d[i][j]*math.asin(norm_stable[i][j]),3) 

# 选取离群度前十名
all_number = []
for item in result:
    for x in item:
        all_number.append(x)
sortedList = sorted(all_number)
result = sortedList[50:]
print(result)
