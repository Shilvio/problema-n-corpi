# settare timout esecuzione
import subprocess as sp
numberBody = 1
for i in range(1):

    numberBody *= 10
    # Generazione Particelle
    sp.run(['../Generate/particleRand', str(numberBody)])

    # ST chianmo script exaustive single thread
    sp.run(['../Single-Thread/Exaustive/Exaustive'])

    # ST chianmo script Barnes-hut single thread
    sp.run(['../Single-Thread/Barnes-hut/Barnes-hut-Bounding-box'])

    # CUDA chiamo script exaustive
    # sp.call(['../CUDA/Exaustive/ExaustiveArrays'])

    # CUDA chiamo script Barnes-hut
    # sp.call(['../CUDA/Barnes-hut/Barnes-hut-bottom-up'])

    # MPI chiamo script Exaustive
    sp.run(['mpiexec', '-n', '5', '../MPI/Exaustive/Exaustive'])

    # MPI chiamo script Barnes-hut
    sp.run(['mpiexec', '-n', '5', '../MPI/Barnes-hut/Barnes-hut'])
