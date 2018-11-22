#-*- coding: utf-8 -*-
from __future__ import division
import arcpy
import numpy as np
import time
import math

"""
算法总体思路：
1、首先根据voronoi多边形以及选择的N-级邻域确定局部区域，确定好局部区域的同时还需要建立一个空间距离矩阵D
2、获取数据后进行数据的归一化
3、计算目标的属性值与周围邻居的差距【注】此处使用的是直接计算局部区域所有的值而不模仿SLOM的距离定义
4、计算区域的信息熵（区域还要再进行一次比例计算）local_entropy，
    同时，根据空间距离矩阵来计算基于空间对象分布的信息熵，用于表征空间分布的稳定性【注】此处使用的是角度信息熵不是距离矩阵
5、最后计算目标的异常值 degree = d-cut/local_entropy 
"""


# 设置工作空间
path = "G:/ArcGIS_Data/privateDataBase.mdb"
arcpy.env.workspace = arcpy.env.workspace = path

"""
定义一个数据表用于存储所有的地理数据，每一条都是以字典的形式存储，格式为：
{
    "RECORD_ID":xx,
    "LM88_POPM88":xxx,
    "Polygon":xxxx,
    "XY":[x,y],       要素的质心xy坐标,用于构建局部区域矩阵
    "norm_data":标准化LM88_POPM88
    "neighbours":[{"RECORD_ID":yy,"XY":[x,y],"LM88_POPM88":yyy,"norm_data":yyyy},{"RECORD_ID":yy2,"LM88_POPM88":yyy2,"norm_data":yyyy2},...]  
}
"""
# 创建一个表格用于存储所有的数据，每一项数据都是要素对应的字典
globe_table = []

with arcpy.da.SearchCursor("美国肺癌数据",["RECORD_ID","LM88_POPM88","SHAPE@","SHAPE@XY"]) as cursor:
    for row in sorted(cursor):
        temp = {}
        temp["RECORD_ID"] = row[0]
        temp["LM88_POPM88"] = row[1]
        temp["Polygen"] = row[2]
        temp["XY"] = row[3]
        # 将要素对象格式化为字典后存储到自己全表中
        globe_table.append(temp)

# 进行数据归一化操作，并且将每个要素归一化的结果存储值字典的“norm_data”中
attrData = []
for item in globe_table:
    attrData.append(item["LM88_POPM88"])
GlobeMax = max(attrData)  # 指定属性中全局最大值
GlobeMin = min(attrData)  # 指定属性中全局最小值
for item in globe_table:
    item["norm_data"] = (item["LM88_POPM88"] - GlobeMin)/(GlobeMax - GlobeMin)
    # 当归一化操作遇到最小值的时候赋一个非常小的数值给它，防止出现除0错误的出现
    if(item["norm_data"]==0):
        item["norm_data"] = (GlobeMax - GlobeMin)*(1/10000000000)

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
            modify["RECORD_ID"] = t["RECORD_ID"] 
            modify["LM88_POPM88"] = t["LM88_POPM88"] 
            modify["norm_data"] = t["norm_data"]
            modify["XY"] = t["XY"]
            neighbours.append(modify)
    return neighbours

# 找到全部对象的一级邻居
for t in globe_table:
    neighb = getNeighbours(t, globe_table)
    t["neighbours"] = neighb


def getAngle(centerPoint, refPoint):
    xi = refPoint[0]
    xo = centerPoint[0]
    yi = refPoint[1]
    yo = centerPoint[1]
    theta = None
    if(xi > xo):
        if(yi > yo):
            theta = math.atan(abs((yi-yo)/(xi-xo)))
        else:
            if(yi == 0):
                theta = 0
            else:
                theta = 2*math.pi - math.atan(abs((yi-yo)/(xi-xo)))
    else:
        if(xi == 0):
            if(yi>yo):
                theta = math.pi/2
            else:
                theta = 3*math.pi/2
        else:
            if(yi>yo):
                theta = math.pi - math.atan(abs((yi-yo)/(xi-xo)))
            else:
                if(yi>yo):
                    theta = math.pi
                else:
                    theta = math.pi + math.atan(abs((yi-yo)/(xi-xo)))
    return theta

# 已知一组点，第一个是参考点，其余的是需要计算角度的点。返回各自的夹角，分摊360
def point2angle(target):
    pointset = [target["XY"]]
    for item in target["neighbours"]:
        pointset.append(item["XY"])
    tempangleSet = []
    for i in range(1,len(pointset)):
        tempangleSet.append(getAngle(pointset[0],pointset[i]))
    # sort all angles
    tempangleSet_sorted = sorted(tempangleSet)
    tempangleSet_sorted.append(2*math.pi+tempangleSet_sorted[0])
    result = []
    for i in range(len(tempangleSet_sorted)-1):
        result.append(tempangleSet_sorted[i+1]-tempangleSet_sorted[i])
    return result

"""
该处使用SLOM中计算属性相对距离作为本文的距离
自己的思路：
1、对局部区域数值进行排序，除去最大值和最小值后计算剩下项的均值
2、计算目标对象和上一步均值的差值作为这个属性相对差距
"""
def getAbsDistance(feature):
    feature_neighbers = feature["neighbours"]
    attr_array = [feature["norm_data"]]
    for item in feature_neighbers:
        attr_array.append(item["norm_data"])
    
    # 除去一个最大值一个最小值
    Max_value = max(attr_array)
    Min_value = min(attr_array)
    attr_array.remove(Max_value)
    attr_array.remove(Min_value)
    avg = sum(attr_array)/len(attr_array)
    # 定义对象的局部距离为数值与周围均值的差异
    d = abs(feature["norm_data"]-avg)
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

        # print("item:%10.7f,sum:%f,attr_array:%d"%(item,sumValue,len(attr_array)))

        AttributeEntropy += (item/sumValue)*math.log(sumValue/item,2)
    return AttributeEntropy

"""
基于产生角度占比计算局部区域的空间分布的信息熵
"""
def getAngleEntropy(feature):
    angleset = point2angle(feature)
    COUNT = len(angleset)
    sumValue = 2*math.pi
    AngleEnrtopy = 0
    for item in angleset:
        AngleEnrtopy += (item/sumValue)*math.log(sumValue/item,2)
    return AngleEnrtopy

# 按照论文方法最后得到的异常值结果
def getScore(feature):
    attributeEntropy = getAttributeEntropy(feature)
    angleEntropy = getAngleEntropy(feature)
    absDistance = getAbsDistance(feature)
    #result = absDistance*(attributeEntropy*0.5+angleEntropy*0.5)
    result = absDistance*attributeEntropy # 首先不加入空间分布信息熵因子，只考虑绝对差值和属性信息熵

    print("paper_d:%f,ArrtributeEntropy:%f,AngleEnrtopy:%f"%(absDistance,attributeEntropy,angleEntropy))
    return result

for item in globe_table:
    item["outlierScore"] = getScore(item)



#创建图层，并且重新将ID、outlierScore、以及geometry属性重新赋值之后再绘制
now = time.strftime('%y%m%d%H%M%S',time.localtime())  # 建立一个时间戳用于命名图层，图层存在的时候重新创建会失败
resultLayer = arcpy.CreateFeatureclass_management(path,"testPaperAlg_2"+now,'POLYGON')

# 添加属性
resultLayer = arcpy.AddField_management(resultLayer,"RECORD_ID","SHORT")
resultLayer = arcpy.AddField_management(resultLayer,"norm_data","DOUBLE")
resultLayer = arcpy.AddField_management(resultLayer,"outlierScore","DOUBLE")




all_fields = ["RECORD_ID","norm_data","outlierScore","SHAPE@"]

# 新建集合用于存储所有的计算结果与几何图形，最后用于重新绘制图层
data = []
# 将所有的点转化成Arcpy格式的Point，然后组合成Arcpy.Array,打包赋值给Arcpy.Polygon()对象
for item in globe_table:
    TEMPitem = [item["RECORD_ID"],item["norm_data"],item["outlierScore"],item["Polygen"]]
    data.append(TEMPitem) 

with arcpy.da.InsertCursor(resultLayer,all_fields) as cursor:
    for row in data:
        cursor.insertRow(row)
print "DONE!"