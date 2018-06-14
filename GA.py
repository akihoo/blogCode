#author:微石
#blog:akihoo.github.io
import numpy as np
import matplotlib.pyplot as plt
import math
plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
def aimFunction(x):
    x=x/2**17*9#将染色体映射到解空间
    y=x+5*math.sin(5*x)+2*math.cos(3*x)
    return y

def Fitness(population):
    value=[]
    for i in range(len(population)):
        value.append(max(aimFunction(population[i]),0))
    return value

def selection(population):
    #计算适应度函数
    fitness=Fitness(population)
    #轮盘赌选择
    sum_fitness=sum(fitness)#计算总体适应度
    fitness=[i/sum_fitness for i in fitness]#计算个体生存概率
    for i in range(1,len(fitness)):
        fitness[i]+=fitness[i-1]
    #优胜劣汰，重新选择10个种群，生存概率高的可能被选择多次
    population_new=[]
    for i in range(len(fitness)):
        rand=np.random.uniform(0,1)
        for j in range(len(fitness)):
            if rand<=fitness[j]:
                population_new.append(population[j]) 
                break          
    return population_new

def get_last_chrom(parent,copint):#计算染色体后几位
    s=''
    for i in range(copint):
        if(parent%2):
            s+='1'
        else:
            s+='0'
        parent=parent>>1
    return int(s,2)

def crossover(population_new,chrom_length, pc=0.6):
    #pc:交配概率
    #chrom_length:染色体长度
    #选择父辈、母辈，交叉染色体
    half=len(population_new)//2
    father,mother=population_new[:half],population_new[half:]
    np.random.shuffle(father)
    np.random.shuffle(mother)
    offspring=[]
    for i in range(half):
        if np.random.uniform(0,1)<=pc:
            copint = np.random.randint(1,chrom_length-1)
            son=(father[i]>>copint<<copint)+get_last_chrom(mother[i],copint)
            daughter=(mother[i]>>copint<<copint)+get_last_chrom(mother[i],copint)
        else:
            son,daughter=father[i],mother[i]
        offspring.append(son)
        offspring.append(daughter)
    return offspring

def mutation(offspring,chrom_length,pm=0.01):
    #pm:变异概率
    for i in range(len(offspring)):
        if np.random.uniform(0,1)<=pm:
            position=np.random.randint(0,chrom_length-1)
            offspring[i]=list(bin(offspring[i])[2:])
            offspring[i]=['0' * (chrom_length-len(offspring))]+offspring[i]
            offspring[i][position]='1' if offspring[i][position]=='0' else '0'
            offspring[i]=int(''.join(offspring[i]),2)
    return offspring

def one_generation(population,chrom_length):
    population=selection(population)#选择
    offspring=crossover(population,chrom_length)#交叉
    offspring=mutation(offspring,chrom_length)#变异
    return max([aimFunction(i) for i in offspring]),offspring

if __name__=='__main__':
    #初始化种群，生成100个种群，使用np.random.randint随机生成基因
    population=[] #储存基因数据
    pop_size=100 #种群规模
    chrom_length=17#基因长度
    for i  in range(pop_size):#随机生成种群
        population.append(np.random.randint(0,2**chrom_length))
    res,pop=[],[]
    for i in range(100):
        r,population=one_generation(population,chrom_length)
        res.append(r)
        if(i%5==0):
            pop.append(population)
    plt.plot(res) 
    plt.show()        
    x=[i/100 for i in range(900)]
    y=[i+5*math.sin(5*i)+2*math.cos(3*i) for i in x]
    colors=['b','r','g','y']
    for j in range(4):
        plt.subplot(221+j)
        plt.plot(x,y)
        plt.scatter([i/2**17*9 for i in pop[j]], [aimFunction(i) for i in pop[j]],color=colors[j],alpha=0.5)
        plt.title('第'+str(j*5)+'代种群')#显示图表标题
    plt.show()
