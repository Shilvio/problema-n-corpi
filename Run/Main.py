# settare timout esecuzione
import subprocess as sp
import resource
import random
numberBody = 1
times = [[], [], [], [], [], []]
for i in range(1):

    numberBody *= 10

    for j in range(1):

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
        # stTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        # sp.run(['../Single-Thread/Barnes-hut/Barnes-hut-Bounding-box'])
        # stTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        # stBarnesHutTime = stTimerEndBH.ru_utime - stTimerStartBH.ru_utime

        # times[1].append(stBarnesHutTime)

        # CUDA chiamo script exaustive
        # cudaTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
        # sp.call(['../CUDA/Exaustive/ExaustiveArrays'])
        # cudaTimerEndE = resource.getrusage(resource.RUSAGE_CHILDREN)
        # cudaExaustiveTime = cudaTimerEnd.ru_utime - cudaTimerStartE.ru_utime

        # times[2].append(cudaExaustiveTime)

        # CUDA chiamo script Barnes-hut
        # cudaTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        # sp.call(['../CUDA/Barnes-hut/Barnes-hut-bottom-up'])
        # cudaTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        # cudaBarnesHutTime = cudaTimerEndBH.ru_utime - cudaTimerStartBH.ru_utime

        # times[3].append(cudaBarnesHutTime)

        # MPI chiamo script Exaustive
        mpiTimerStartE = resource.getrusage(resource.RUSAGE_CHILDREN)
        sp.run(['mpiexec', '-n', '5', '../MPI/Exaustive/Exaustive'])
        mpiTimerEndE = resource.getrusage(resource.RUSAGE_CHILDREN)
        mpiExaustiveTime = mpiTimerEndE.ru_utime - mpiTimerStartE.ru_utime

        times[4].append(mpiExaustiveTime)

        # MPI chiamo script Barnes-hut
        # mpiTimerStartBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        # sp.run(['mpiexec', '-n', '5', '../MPI/Barnes-hut/Barnes-hut'])
        # mpiTimerEndBH = resource.getrusage(resource.RUSAGE_CHILDREN)
        # mpiBarnesHutTime = mpiTimerEndBH.ru_utime - mpiTimerStartBH.ru_utime

        # times[5].append(mpiBarnesHutTime)

    print(times[0])
    print(times[4])
