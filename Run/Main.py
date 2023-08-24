# settare timout esecuzione
import subprocess as sp
import resource
import random
#staring variables
numberBody = 1
timesBody = [[], [], [], [], [], []]
timeIter=[]

# iteration, mean and subprocess timeout
bodyIterations = 1
meanIterations = 5
timeout_s = 20

print()
for i in range(bodyIterations):

    numberBody *= 10

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
                times[0].append(-1)
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
                times[1].append(-1)
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
                times[2].append(-1)
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
                times[3].append(-1)
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
                times[4].append(-1)
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
                times[5].append(-1)
                print(f'Timeout for mpi Barnes-hut ({timeout_s}s) expired')
                mpiB=False

    print()
    timeIter.append(times)

for bodys in timeIter:
    count=0
    for time in bodys:
        temp=0
        for exec in time:
            if(exec==-1):
                temp= (-len(time))
                break
            else:
                temp+=exec
        if len(time)!=0:
            timesBody[count].append(temp/len(time))
        count+=1

print("st exaustive:")
print(timesBody[0])
print("st Barnes-hut:")
print(timesBody[1])
print()
print("Cuda exaustive:")
print(timesBody[2])
print("Cuda Barnes-hut:")
print(timesBody[3])
print()
print("MPI exaustive:")
print(timesBody[4])
print("MPI Barnes-hut:")
print(timesBody[5])
