# settare timout esecuzione
import subprocess as sp
numberBody = 1
for i in range(2):

    numberBody *= 10
    sp.call(['../Generate/particleRand', str(numberBody)])

    # CHIAMO SCRIPT SINGLE THREAD
    sp.call(['../Single-Thread/Exaustive/Exaustive'])

    # CHIAMO SCRIPT SINGLE THREAD
    sp.call(['../Single-Thread/Barnes-hut/Barnes-hut'])
