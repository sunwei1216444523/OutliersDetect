#-*- coding: utf-8 -*-
from __future__ import division
import arcpy
import numpy as np
import time


# 设置工作空间
arcpy.env.workspace = "G:/ArcGIS_Data/privateDataBase.mdb"
path = "G:/ArcGIS_Data/privateDataBase.mdb"
"""
定义一个数据表用于存储所有的地理数据，每一条都是以字典的形式存储，格式为：
{
    "RECORD_ID":xx,
    "LM88":xxx,
    "Polygen":xxxx,
    "norm_data":标准化LM88
    "neighbours":[{"RECORD_ID":yy,"LM88":yyy,"norm_data":yyyy},{"RECORD_ID":yy2,"LM88":yyy2,"norm_data":yyyy2},...]  
}
"""
all_table = []
# 测试数据为：美国肺癌数据
# 使用with语句创建游标，注意此处文件是没有后缀的，加后缀读取不了
with arcpy.da.SearchCursor("美国肺癌数据", ["RECORD_ID","LM88","SHAPE@"]) as cursor:
    for row in sorted(cursor):
        temp = {}
        temp["RECORD_ID"] = row[0]
        temp["LM88"] = row[1]
        temp["Polygen"] = row[2]
        # 将所有表格数据保存至字典中，随后字典中还需要 保存每个对象的邻居
        all_table.append(temp)


# 数据的归一化应该放在全局中执行，而不是局部归一化
mydata = []
for item in all_table:
    mydata.append(item["LM88"])
mmax = max(mydata)
mmin = min(mydata)
print '-----------------------------'
print("max:",mmax)
print("min:",mmin)
print '--------------------------------'
for item in all_table:
    item["norm_data"] = (item["LM88"]-mmin)/(mmax-mmin)

#===================================================
"""
该函数用于获取某个对象的所有邻居
target 目标对象
all_data 图中所有对象
return 目标对邻居的集合
"""
def getNeighbours(target, all_data):
    neighbours = []
    for t in all_data:
        if t["Polygen"].touches(target["Polygen"]):
            # 只保存邻居的ID 和 要计算的数据，后期有其他需求再添加
            modify = {}
            modify["RECORD_ID"] = t["RECORD_ID"] 
            modify["LM88"] = t["LM88"] 
            modify["norm_data"] = t["norm_data"]
            neighbours.append(modify)
    return neighbours

# 找到所有空间实体的对象
for t in all_table:
    neighb = getNeighbours(t, all_table)
    t["neighbours"] = neighb

# ======测试一号对象的邻居有：2,7,6,11============
# one = all_table[0]
# n = one["neighbours"] # n是一个字典的集合
# for i in n :
#     print(i["RECORD_ID"])
#================================================


# 每个对象的邻居都保存在各自的数据中，接下来就是算法的实现

# 算法一：SLOM
"""
1、数据的归一化，策略是将所有的数据放在一个数组，目标数据放在第一个
2、
"""
# def normalization(allValue):
#     maxm = np.max(allValue)
#     minm = np.min(allValue)
#     result = []
#     for i in allValue:
#         r = (i-minm)/(maxm-minm)
#         result.append(r)
#     return result

"""
计算SLOM算法中定义的d值
"""
def getSLOM_d(o):
    o_neighbs = o["neighbours"]
    # 将局部所有的数据都放入一个数组中,目标对象放在第一个
    attr = [o["norm_data"]]
    for ii in o_neighbs:
        attr.append(ii["norm_data"])
    # 数据归一化
    # norm_data = normalization(attr)

    # N(o)
    N_o = len(attr) - 1
    # 所有邻居到目标的差距
    diff  = []
    for i in range(len(attr)):
        if i == 0:
            continue 
        diff.append(abs(attr[i]-attr[0]))
    if o["RECORD_ID"]==77:
        print "-=-=-=-=-=test demo-=-=-=-=-=-=-=-"
        print attr
        print diff
        print "-=-=-=-=-=-=-=--=-==-=-=-=-=-"

    total_diff = sum(diff)
    max_diff = max(diff)
    
    # 计算SLOM_d
    SLOM_d = (total_diff-max_diff)/(N_o-1)
    return SLOM_d

"""
计算局部的SLOM中的d值
p 表示衷心计算的点
all_local_object 表示该区域下所有的目标
"""
def getSLOM_d_local(p,all_loacal_objct):
    SLOM_d_loacl = 0
    # 获取点p所有的邻居
    p_neighbours = p["neighbours"]
    
    # 先定义一个局部邻居的缓存集合
    local_tamp = []
    # 得到区域下p的所有邻居，用全局邻居和局部进行求交集
    for pp in p_neighbours:
        for ppp in all_loacal_objct:
            if pp["RECORD_ID"] == ppp["RECORD_ID"]:
                local_tamp.append(pp)

    # test33 = []
    # for test3 in local_tamp:
    #     test33.append(test3["RECORD_ID"])
    # print("subset point include",test33)

    # 开始计算实际数值
    center_num = p["norm_data"] 
    conv2num = []
    for xx in local_tamp:
        conv2num.append(xx["norm_data"])
    
    diff_array = []
    for i in conv2num:
        diff_array.append(abs(i-center_num))
    
    SLOM_d_loacl = (sum(diff_array)-max(diff_array))/(len(local_tamp)+1)
    return SLOM_d_loacl

"""
计算震荡因子β
"""
def getSLOM_beita(o):
    # 初始化β
    SLOM_beita = 0

    #  区域内所有的对象,包括中心点
    tempObject = [o]

    o_neighbs = o["neighbours"]

    # 将局部所有归一化之后的数据都放入一个数组中，包括中心点
    attr = [o["norm_data"]]

    # 因为每个对象的“neighbours”中的邻居信息只有ID和属性，所以需要根据ID找到完整的邻居信息
    for ii in o_neighbs:
        attr.append(ii["norm_data"])
        for iii in all_table:
            if ii["RECORD_ID"]==iii["RECORD_ID"]:
                tempObject.append(iii)

    #N+(o)
    N_plus_o = len(tempObject)
    s = 0
    for p in tempObject:
        t = getSLOM_d_local(p,tempObject)
        s += t
    
    # avg(N+(o))
    avg_N_plus_o = s/N_plus_o

    for p in tempObject:
        if getSLOM_d_local(p,tempObject) > avg_N_plus_o:
            SLOM_beita += 1
        else: 
            if getSLOM_d_local(p,tempObject) < avg_N_plus_o:
                SLOM_beita -= 1

    SLOM_beita  = abs(SLOM_beita)
    # 该部分的计算是SLOM 的一个限制，当邻居只有一个的时候，就会出现分母为0的情况，所以加上一个极小的数
    SLOM_beita = max([SLOM_beita,1])/(N_plus_o-2+0.000000001)
    SLOM_beita = SLOM_beita/(1+(s/(len(attr)-1)))

    result = SLOM_beita
    return result

def get_SLOM(o):
    print(o["RECORD_ID"],"d:",getSLOM_d(o),"beita:",getSLOM_beita(o),"SLOM:",getSLOM_d(o)*getSLOM_beita(o))
    return getSLOM_d(o)*getSLOM_beita(o)


  
# 为每个对象计算SLOM值
for item in all_table:
    item["SLOM"] = get_SLOM(item) or 0;


print "============================================"
# 打印所有的数据的归一化数据以及SLOM数据
for x in all_table:
    print("norm_data:",x["norm_data"],"SLOM:",x["SLOM"])
#创建图层，并且重新将ID、SLOM、以及geometry属性重新赋值之后再绘制
now = time.strftime('%y%m%d%H%M%S',time.localtime())  # 建立一个时间戳用于命名图层，图层存在的时候重新创建会失败
resultLayer = arcpy.CreateFeatureclass_management(path,"testSLOM"+now,'POLYGON')

# 添加属性
resultLayer = arcpy.AddField_management(resultLayer,"RECORD_ID","SHORT")
resultLayer = arcpy.AddField_management(resultLayer,"norm_data","DOUBLE")
resultLayer = arcpy.AddField_management(resultLayer,"SLOM","DOUBLE")


all_fields = ["RECORD_ID","norm_data","SLOM","SHAPE@"]

# 新建集合用于存储所有的计算结果与几何图形，最后用于重新绘制图层
data = []
# 将所有的点转化成Arcpy格式的Point，然后组合成Arcpy.Array,打包赋值给Arcpy.Polygon()对象
for item in all_table:
    TEMPitem = [item["RECORD_ID"],item["norm_data"],item["SLOM"],item["Polygen"]]
    data.append(TEMPitem) 

with arcpy.da.InsertCursor(resultLayer,all_fields) as cursor:
    for row in data:
        cursor.insertRow(row)
print "DONE!"







    
    

    

