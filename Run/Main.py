# settare timout esecuzione
import subprocess as sp
import resource
import random
import matplotlib.pyplot as plt
#staring variables
timesBody = [[], [], [], [], [], []]
timeIter=[]
timeout_s = 60

# tempi di esecuzione degli algoritmi
def timing():
    
    numberBody=1
    nbody=[]
    bodyIterations = 6
    meanIterations = 5

    print()
    for i in range(bodyIterations):

        numberBody *= 10
        nbody.append(numberBody)

        stE=True;stB=True
        cudaE=True;cudaB=True
        mpiE=True;mpiB=True

        print("\n"+str(numberBody) +"\n")
        times = [[], [], [], [], [], []]
        for j in range(meanIterations):
            
            print(j+1)

            seed = random.randint(-5*10E3, 5*10E3)

            # Generazione Particelle
            sp.run(['../Generate/particleRand', str(numberBody), str(seed)])

            # ST chianmo script exaustive single thread
            if stE:
                try:
                    stTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
                    sp.run(['../Single-Thread/Exaustive/Exaustive'],timeout=timeout_s)
                    stTimerEndE = resource.getrusage(resource.RUSAGE_CHILDREN)
                    stExaustiveTime = stTimerEndE.ru_utime - stTimerStartE.ru_utime

                    times[0].append(stExaustiveTime)
                    
                except sp.TimeoutExpired:
                    times[0].append(timeout_s+1)
                    print(f'Timeout for st exaustive ({timeout_s}s) expired')
                    stE=False
                
            # ST chianmo script Barnes-hut single thread
            if stB:
                try:
                    stTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
                    sp.run(['../Single-Thread/Barnes-Hut/Barnes-hut-Bounding-box'],timeout=timeout_s)
                    stTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
                    stBarnesHutTime = stTimerEndBH.ru_utime - stTimerStartBH.ru_utime

                    times[1].append(stBarnesHutTime)

                except sp.TimeoutExpired:
                    times[1].append(timeout_s+1)
                    print(f'Timeout for st Barnes-hut ({timeout_s}s) expired')
                    stB=False

            # CUDA chiamo script exaustive
            if cudaE:
                try:
                    cudaTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
                    sp.call(['../CUDA/Exaustive/ExaustiveArrays'],timeout=timeout_s)
                    cudaTimerEnd = resource.getrusage(resource.RUSAGE_CHILDREN)
                    cudaExaustiveTime = cudaTimerEnd.ru_utime - cudaTimerStartE.ru_utime

                    times[2].append(cudaExaustiveTime)
                    
                except sp.TimeoutExpired:
                    times[2].append(timeout_s+1)
                    print(f'Timeout for cuda exaustive ({timeout_s}s) expired')
                    cudaE=False

            # CUDA chiamo script Barnes-hut
            if cudaB:
                try:
                    cudaTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
                    sp.call(['../CUDA/Barnes-Hut/Barnes-hut-bottom-up'],timeout=timeout_s)
                    cudaTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
                    cudaBarnesHutTime = cudaTimerEndBH.ru_utime - cudaTimerStartBH.ru_utime  
                    
                    times[3].append(cudaBarnesHutTime)
                    
                except sp.TimeoutExpired:
                    times[3].append(timeout_s+1)
                    print(f'Timeout for cuda Barnes-hut ({timeout_s}s) expired')
                    cudaB=False

            # MPI chiamo script Exaustive
            if mpiE:
                try:
                    mpiTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
                    sp.run(['mpiexec', '-n', '6', '../MPI/Exaustive/Exaustive'],timeout=timeout_s)
                    mpiTimerEndE = resource.getrusage(resource.RUSAGE_CHILDREN)
                    mpiExaustiveTime = mpiTimerEndE.ru_utime - mpiTimerStartE.ru_utime

                    times[4].append(mpiExaustiveTime)
                    
                except sp.TimeoutExpired:
                    times[4].append(timeout_s+1)
                    print(f'Timeout for mpi exaustive ({timeout_s}s) expired')
                    mpiE=False
            
            # MPI chiamo script Barnes-hut
            if mpiB:
                try:
                    mpiTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
                    sp.run(['mpiexec', '-n', '6', '../MPI/Barnes-Hut/Barnes-hut'],timeout=timeout_s)
                    mpiTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
                    mpiBarnesHutTime = mpiTimerEndBH.ru_utime - mpiTimerStartBH.ru_utime

                    times[5].append(mpiBarnesHutTime)
                    
                except sp.TimeoutExpired:
                    times[5].append(timeout_s+1)
                    print(f'Timeout for mpi Barnes-hut ({timeout_s}s) expired')
                    mpiB=False

        print()
        timeIter.append(times)

    for bodys in timeIter:
        count=0
        for time in bodys:
            temp=0
            for exec in time:
                temp+=exec
            if len(time)!=0:
                timesBody[count].append(temp/len(time))
            count+=1

    print(timesBody)
    print(nbody)
    return timesBody,nbody

def graphGenerator(timeData,nbody):

    plt.xlim(left=0)
    plt.xlim(right=timeout_s)
    
    #stE
    plt.plot(timeData[0],nbody,'ro-',label='st exaustive')
    #stB
    plt.plot(timeData[1],nbody,'ro--',label='st barnes-hut')
    #cudaE
    plt.plot(timeData[2],nbody,'go-',label='cuda exaustive')
    #cudaB
    plt.plot(timeData[3],nbody,'go--',label='cuda barnes-hut')
    #mpiE
    plt.plot(timeData[4],nbody,'bo-',label='mpi exaustive')
    #mpiB
    plt.plot(timeData[5],nbody,'bo--',label='mpi barnes-hut')

    plt.yscale('log')
    plt.xscale('function', functions=(lambda x: x**0.3, lambda x: x**3))
    plt.grid(True)
    plt.ylabel('Body')
    plt.xlabel('Time')
    plt.legend(loc='lower right')
    
    plt.show()
    
    pass

if __name__ == '__main__':
    #timeData,nbody=timing()

    tempData = [[0.0005625999999999813, 0.0008181999999999912, 0.07054760000000009, 7.1653272, 61.0, 61.0], [0.00046640000000000933, 0.0015291999999999862, 0.010622199999999981, 0.17193819999999818, 3.153712799999994, 50.30531179999999], [0.005809999999999999, 0.00799840000000005, 0.035531200000000054, 0.3253437999999971, 8.821345000000008, 61.0], [0.008587200000000008, 0.012182399999999972, 0.01841779999999984, 0.10229660000000287, 0.46098200000000134, 4.171039000000019], [0.067489, 0.0617316, 0.18768840000000003, 11.182502000000003, 61.0, 61.0], [0.0720584, 0.06659540000000001, 0.11230300000000008, 0.7754125999999999, 10.72308099999999, 133.8670334]]
    tempBody= [10, 100, 1000, 10000, 100000, 1000000]


    graphGenerator(tempData,tempBody)

