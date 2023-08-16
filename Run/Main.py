# settare timout esecuzione
import subprocess as sp
numberBody = 1
for i in range(2):

    numberBody *= 10
    sp.call(['../Generate/particleRand', str(numberBody)])

    

    # ST chianmo script exaustive single thread
    sp.call(['../Single-Thread/Exaustive/Exaustive'])

    # ST chianmo script Barnes-hut single thread
    sp.call(['../Single-Thread/Barnes-hut/Barnes-hut'])


    # CUDA chiamo script exaustive 
    sp.call(['../CUDA/Exaustive/ExaustiveArrays'])

    # CUDA chiamo script Barnes-hut 
    sp.call(['../CUDA/Barnes-hut/Barnes-hut-bottom-up'])


    # MPI chiamo script Exaustive 
    #sp.call(['../MPI/Exaustive/Exaistive'])

    # MPI chiamo script Barnes-hut 
    #sp.call(['../MPI/Barnes-hut'])
