#-*- coding: utf-8 -*-
from __future__ import division
import arcpy
import numpy as np
import time
import math


# 设置工作空间
path = "G:/ArcGIS_Data/privateDataBase.mdb"
arcpy.env.workspace = arcpy.env.workspace = path
"""
定义一个数据表用于存储所有的地理数据，每一条都是以字典的形式存储，格式为：
{
    "NAME":xx,
    "RECORD_ID":xx,
    "LM88_POPM88":xxx,
    "Polygon":xxxx,
    "XY":[x,y],       要素的质心xy坐标,用于构建局部区域矩阵
    "norm_data":标准化LM88_POPM88
    "neighbours":[{"NAME":y,RECORD_ID":yy,"XY":[x,y],"LM88_POPM88":yyy,"norm_data":yyyy},{"NAME":y,"RECORD_ID":yy2,"LM88_POPM88":yyy2,"norm_data":yyyy2},...]  
}
"""

"""
该函数用于获取某个对象的所有邻居
target 目标对象
globe_table 图中所有对象
return 目标对邻居的集合
"""
def getNeighbours(target, globe_table):
    neighbours = []
    for t in globe_table:
        if t["Polygen"].touches(target["Polygen"]):
            # 只保存邻居的ID 和 要计算的数据，后期有其他需求再添加
            modify = {}
            modify["NAME"] = t["NAME"]
            modify["RECORD_ID"] = t["RECORD_ID"] 
            modify["LM88_POPM88"] = t["LM88_POPM88"] 
            modify["norm_data"] = t["norm_data"]
            modify["XY"] = t["XY"]
            neighbours.append(modify)
    return neighbours

"""
update:
局部区域只去除差异值最大的目标，类似于SLOM的定义
"""
def getAbsDistance(feature):
    feature_neighbers = feature["neighbours"]
    diff_array = []
    for item in feature_neighbers:
        diff_array.append(abs(feature["norm_data"]-item["norm_data"]))
    # 定义对象的局部距离为数值与周围均值的差异
    d = (sum(diff_array)-max(diff_array))/(len(diff_array)-1)
    return d

"""
计算局部要素属性值的信息熵
"""
def getAttributeEntropy(feature):
    feature_neighbers = feature["neighbours"]
    attr_array = []
    sumValue = 0
    for item in feature_neighbers:
        attr_array.append(item["norm_data"])
        sumValue += item["norm_data"]
    # 以属性值占比作为概率,计算信息熵
    AttributeEntropy = 0
    for item in attr_array:
        p = (item+0.0000001)/(sumValue+0.0000001*len(attr_array))
        AttributeEntropy += p*math.log(1/p,2)
    rate_entropy = AttributeEntropy/math.log(len(feature_neighbers),2)
    return rate_entropy

# 按照论文方法最后得到的异常值结果
def getScore(feature):
    result = feature["d"]*math.asin(feature["norm_stable"])
    return result

# 创建一个表格用于存储所有的数据，每一项数据都是要素对应的字典
globe_table = []

with arcpy.da.SearchCursor("美国肺癌数据",["NAME","RECORD_ID","LM88_POPM88","SHAPE@","SHAPE@XY"]) as cursor:
    for row in sorted(cursor):
        temp = {}
        temp["NAME"] = row[0]
        temp["RECORD_ID"] = row[1]
        temp["LM88_POPM88"] = row[2]
        temp["Polygen"] = row[3]
        temp["XY"] = row[4]
        # 将要素对象格式化为字典后存储到自己全表中
        globe_table.append(temp)

# 进行数据归一化操作，并且将每个要素归一化的结果存储值字典的“norm_data”中
attrData = []
#  源数据归一化后结果的一维数组
DataOneDim = [] 
for item in globe_table:
    attrData.append(item["LM88_POPM88"])

# 指定属性中全局最大值
GlobeMax = max(attrData)  
# 指定属性中全局最小值  
GlobeMin = min(attrData) 


for item in globe_table:
    item["norm_data"] = (item["LM88_POPM88"] - GlobeMin)/(GlobeMax - GlobeMin)
    DataOneDim.append(item["norm_data"])



# 找到全部对象的一级邻居
for t in globe_table:
    neighb = getNeighbours(t, globe_table)
    t["neighbours"] = neighb

# 计算出稳定性之后，还需要归一化
all_stable_value = []
for item in globe_table:
    item["stable"] = getAttributeEntropy(item)
    all_stable_value.append(item["stable"])
max_stable = max(all_stable_value)
min_stable = min(all_stable_value)
for item in globe_table:
    # 归一化把稳定度拉伸至[0.1-1]之间
    item["norm_stable"] = (item["stable"]-min_stable)/(max_stable-min_stable)*(1-0.1)+0.1

for item in globe_table:
    item["d"] = getAbsDistance(item)

for item in globe_table:
    item["SLOI"] = getScore(item)

#创建图层，并且重新将ID、SLOI、以及geometry属性重新赋值之后再绘制
now = time.strftime('%y%m%d%H%M%S',time.localtime())  # 建立一个时间戳用于命名图层，图层存在的时候重新创建会失败
resultLayer = arcpy.CreateFeatureclass_management(path,"SLOI_arcsin"+now,'POLYGON')

# 添加属性
resultLayer = arcpy.AddField_management(resultLayer,"NAME","TEXT")
resultLayer = arcpy.AddField_management(resultLayer,"RECORD_ID","SHORT")
resultLayer = arcpy.AddField_management(resultLayer,"norm_data","DOUBLE")
resultLayer = arcpy.AddField_management(resultLayer,"d","DOUBLE")
resultLayer = arcpy.AddField_management(resultLayer,"stable","DOUBLE")
resultLayer = arcpy.AddField_management(resultLayer,"norm_stable","DOUBLE")
resultLayer = arcpy.AddField_management(resultLayer,"SLOI","DOUBLE")

all_fields = ["NAME","RECORD_ID","norm_data","d","stable","norm_stable","SLOI","SHAPE@"]

# 新建集合用于存储所有的计算结果与几何图形，最后用于重新绘制图层
data = []
# 将所有的点转化成Arcpy格式的Point，然后组合成Arcpy.Array,打包赋值给Arcpy.Polygon()对象
for item in globe_table:
    TEMPitem = [item["NAME"],item["RECORD_ID"],item["norm_data"],item["d"],item["stable"],item["norm_stable"],item["SLOI"],item["Polygen"]]
    data.append(TEMPitem) 

#写入文件，输出
with arcpy.da.InsertCursor(resultLayer,all_fields) as cursor:
    for row in data:
        cursor.insertRow(row)
print "DONE!"