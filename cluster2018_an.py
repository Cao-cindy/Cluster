#初始数据取平均，再换算成m

import math    
import numpy as np    
import pandas as pd

KN=14 #簇数量

def Exy(n,a):
    #n为聚类的类别数量，a为ɛ二维矩阵
    rows = n
    cols = n
    e = [[0 for j in range(cols)] for i in range(rows)]
    total_sum = 0
    for i in range(0,n):
        for j in range(0,n):
            total_sum = 0
            for k in range(0,n):
                axy = (a[i][k]-a[j][k])*(a[i][k]-a[j][k])
                total_sum += axy
            e[i][j] = total_sum/(4*n)
    return e


def S(positive_m,negative_m):
    #positive_m和negative_m是两个数
    #positive_m是类x与类y之间协同关系的数量，negative_m是类x与类y之间拮抗关系的数量
    if (positive_m*negative_m)!=0:
        positive_p = positive_m/(positive_m+negative_m)
        negative_p= negative_m/(positive_m+negative_m)
        s=(positive_m+negative_m)*( positive_p*math.log(positive_p)+ negative_p*math.log( negative_p))
    elif positive_m == 0 and negative_m!=0:
        negative_p= negative_m/(positive_m+negative_m)
        s=(negative_m)*(negative_p*math.log( negative_p))
    elif negative_m == 0 and positive_m!=0:
        positive_p = positive_m/(positive_m+negative_m)
        s=(positive_m)*( positive_p*math.log(positive_p))
    else:
        s=0
    return s

def Sxy(positive_m,negative_m,n):
    #返回二维数组▲Sxy
    #positive_m和negative_m是两个二维数组,类x和类y之间的协同关系和拮抗关系的数量
    #已知聚类1，2，3，4,5和他们之间的协同和拮抗关系的条数
    rows = n
    cols = n
    s = [[0 for j in range(cols)] for i in range(rows)]
    total_sum = 0
    #计算▲S[i][j]
    for i in range(0,n):
        for j in range(i,n):
            total_sum = 0
            for k in range(0,n):
                if k!=i and k!=j :
                    a = S(positive_m[i][k]+positive_m[j][k],negative_m[i][k]+negative_m[j][k])
                    b = S(positive_m[i][k],negative_m[i][k])
                    c = S(positive_m[j][k],negative_m[j][k])
                    total_sum += abs(a - b - c)
            s[i][j]=S(positive_m[i][j],negative_m[i][j]) - total_sum
            s[j][i] = s[i][j]  # 同时计算 s[i][j] 和 s[j][i]
    return s

def Fxy(antibiotic_dict,cluster1,cluster2,e,s):
    #cluster1,cluster2为cluster[i],cluster[j]，即类x和类y
    #s为▲Sxy
    #e为exy二维数组

    #t为参数，可改
    t=0.1
    mine= float('inf')
    for item1 in cluster1:
        i = antibiotic_dict[item1]
        for item2 in cluster2:
            j = antibiotic_dict[item2]
            mine = min(mine, e[i][j])
    fxy=mine - t*s
    return fxy

def remove_jth_row_col(matrix, j):
    # 删除第j行
    matrix = [row for i, row in enumerate(matrix) if i != j]
    
    # 删除第j列
    matrix = [[row[i] for i in range(len(row)) if i != j] for row in matrix]
    
    return matrix

#纯度
def Purity(positive_m,negative_m):
    #number是簇数量
    sum=0
    sumN=0
    for i in range(KN):
        for j in range(i+1,KN):
            m=max(positive_m[i][j],negative_m[i][j])
            sum=sum+m
            N=positive_m[i][j]+negative_m[i][j]
            sumN=sumN+N
    p=sum/sumN
    return p

#聚类
def Clustering(nodes,exy,antibiotic_dict,positive_m,negative_m,n):
    # 初始化每个元素为一个单独的类别
    clusters = [[node] for node in nodes]
    round=0
    #簇数量
    while len(clusters) > KN: #簇数量
        round=round+1
        print("第{}轮合并".format(round))
        min_distance = float('inf')
        merge_clusters = None
        sxy=Sxy(positive_m,negative_m,len(clusters))
        # 找到距离最近的一对类别
        for i in range(len(clusters)):
            #print("循环到第{}类药物对应行".format(i))
            for j in range(i+1,len(clusters)):

                #if (len(clusters[i]==1) and len(clusters[j]==1)) or sxy[i][j]!=0:
                distance= Fxy(antibiotic_dict,clusters[i],clusters[j],exy,sxy[i][j])
                if distance < min_distance:
                    min_distance = distance
                    #min_mine=mine
                    merge_clusters = (i, j)



        # 合并距离最近的两个类别
        i, j = merge_clusters
        print(f"合并簇{clusters[i]}和{clusters[j]}")
        clusters[i].extend(clusters[j])
        del clusters[j]
        #更新positive_m,negative_m矩阵
        for p in range(0,len(clusters)):
            positive_m[i][p] = positive_m[i][p]+positive_m[j][p]
            positive_m[p][i] = positive_m[i][p]
            negative_m[i][p] = negative_m[i][p]+negative_m[j][p]
            negative_m[p][i] = negative_m[i][p]
        positive_m=remove_jth_row_col(positive_m,j)
        negative_m=remove_jth_row_col(negative_m,j)
    #计算纯度
    p=Purity(positive_m,negative_m) #簇数量
    print("purity:",p)
    return clusters

         

def main():
 # 读取Excel文件
    df = pd.read_excel('D:/研究生/论文/新工作/工作/2018抗生素code.xlsx')

    # 将df转换为二维列表
    nodes = sorted(set(df['Drug1']) | set(df['Drug2']))
    antibiotic_dict = {node: i for i, node in enumerate(nodes)}
    axy = [[None for _ in range(len(nodes))] for _ in range(len(nodes))]
    for _, row in df.iterrows():
        i = nodes.index(row['Drug1'])
        j = nodes.index(row['Drug2'])
        
        e_coli_value = row['PA14']
        if e_coli_value == 'Synergy':
            axy[i][j] = 1
            axy[j][i] = 1
        elif e_coli_value == 'Antagonism':
            axy[i][j] = -1
            axy[j][i] = -1
        else :
            axy[i][j] = 0
            axy[j][i] = 0

    axy = [[0 if val is None else val for val in inner] for inner in axy]

    # 输出更新后的二维列表axy
    
    #print("nodes:",nodes)
    print("antibiotic_dict:",antibiotic_dict)
    print("axy:",axy) 
    

#计算positive_m,negative_m,抗生素之间的初始矩阵
    n = len(axy)
    positive_m = [[0 for _ in range(n)] for _ in range(n)]
    negative_m = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(0,n):
        for j in range(0,n):
            if axy[i][j]  > 0  :
                positive_m[i][j]= 1
                negative_m[i][j]= 0
            elif axy[i][j] < 0 :
                positive_m[i][j]= 0
                negative_m[i][j]= 1
            elif axy[i][j] == 0 :
                positive_m[i][j] = 0
                negative_m[i][j]= 0
    #print("positive_m:",positive_m)  
    #print("negative_m:",negative_m)         

    #计算exy
    exy = [[0 for _ in range(n)] for _ in range(n)]
    exy=Exy(n,axy)
    #print("exy:",exy)  
    '''
    sxy = [[0 for _ in range(n)] for _ in range(n)]
    sxy=Sxy(positive_m,negative_m,n)
    #print("Sxy:",sxy)  
    '''
    clusters=Clustering(nodes,exy,antibiotic_dict,positive_m,negative_m,n)
    print("clusters:",clusters) 

     

if __name__ == "__main__":
    main()