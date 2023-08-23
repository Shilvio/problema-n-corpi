# settare timout esecuzione
import subprocess as sp
import resource
import random
numberBody = 1
times = [[], [], [], [], [], []]
for i in range(3):

    numberBody *= 10

    print("\n"+ str(numberBody) +"\n")
    times = [[], [], [], [], [], []]
    for j in range(5):

        seed = random.randint(-5*10E3, 5*10E3)

        # Generazione Particelle
        sp.run(['../Generate/particleRand', str(numberBody), str(seed)])

        # ST chianmo script exaustive single thread
        stTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
        sp.run(['../Single-Thread/Exaustive/Exaustive'])
        stTimerEndE = resource.getrusage(resource.RUSAGE_CHILDREN)
        stExaustiveTime = stTimerEndE.ru_utime - stTimerStartE.ru_utime

        times[0].append(stExaustiveTime)

        # ST chianmo script Barnes-hut single thread
        stTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        sp.run(['../Single-Thread/Barnes-Hut/Barnes-hut-Bounding-box'])
        stTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        stBarnesHutTime = stTimerEndBH.ru_utime - stTimerStartBH.ru_utime

        times[1].append(stBarnesHutTime)

        # CUDA chiamo script exaustive
        cudaTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
        sp.call(['../CUDA/Exaustive/ExaustiveArrays'])
        cudaTimerEnd = resource.getrusage(resource.RUSAGE_CHILDREN)
        cudaExaustiveTime = cudaTimerEnd.ru_utime - cudaTimerStartE.ru_utime

        times[2].append(cudaExaustiveTime)

        # CUDA chiamo script Barnes-hut
        cudaTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        sp.call(['../CUDA/Barnes-Hut/Barnes-hut-bottom-up'])
        cudaTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        cudaBarnesHutTime = cudaTimerEndBH.ru_utime - cudaTimerStartBH.ru_utime  
        
        times[3].append(cudaBarnesHutTime)

        # MPI chiamo script Exaustive
        mpiTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
        sp.run(['mpiexec', '-n', '6', '../MPI/Exaustive/Exaustive'])
        mpiTimerEndE = resource.getrusage(resource.RUSAGE_CHILDREN)
        mpiExaustiveTime = mpiTimerEndE.ru_utime - mpiTimerStartE.ru_utime

        times[4].append(mpiExaustiveTime)
        


          # MPI chiamo script Barnes-hut
        mpiTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        sp.run(['mpiexec', '-n', '6', '../MPI/Barnes-Hut/Barnes-hut'])
        mpiTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        mpiBarnesHutTime = mpiTimerEndBH.ru_utime - mpiTimerStartBH.ru_utime

        times[5].append(mpiBarnesHutTime)

    print("st exaustive:")
    print(times[0])
    print("st Barnes-hut:")
    print(times[1])
    print("Cuda exaustive:")
    print(times[2])
    print("Cuda Barnes-hut:")
    print(times[3])
    print("MPI exaustive:")
    print(times[4])
    print("MPI Barnes-hut:")
    print(times[5])
